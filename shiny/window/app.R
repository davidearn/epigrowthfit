library("shiny")

ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  titlePanel(
    title = div(HTML("<b>epigrowthfit</b>: Fitting window selection tool")),
    windowTitle = "epigrowthfit: Fitting window selection tool"
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "data",
        label = HTML("Data frame of time series [<tt>.rds</tt> file, required]:"),
        multiple = FALSE,
        accept = ".rds"
      ),
      textInput(
        inputId = "formula",
        label = "Formula indicating data frame structure:",
        value = "cbind(time, x) ~ ts"
      ),
      fileInput(
        inputId = "data_windows",
        label = HTML("Data frame of fitting windows [<tt>.rds</tt> file, optional]:"),
        multiple = FALSE,
        accept = ".rds"
      ),
      textInput(
        inputId = "formula_windows",
        label = "Formula indicating data frame structure:",
        value = "cbind(start, end) ~ ts"
      ),
      tabsetPanel(
        id = "tabset_ts",
        type = "hidden",
        tabPanelBody("tabset_ts_null"),
        tabPanelBody("tabset_ts_select", uiOutput("ui_tabset_ts_select"))
      ),
      uiOutput("ui_tabset_spar")
    ),
    mainPanel(
      fluidRow(
        uiOutput("ui_tabset_plot")
      ),
      fluidRow(
        column(6,
          uiOutput("ui_tabset_points")
        ),
        column(6,
          uiOutput("caption_windows"),
          uiOutput("ui_tabset_windows")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  make_formula <- function(input, id) {
    formula <- try(as.formula(str2lang(input[[id]]), env = .GlobalEnv))
    if (inherits(formula, "try-error")) {
      shinyFeedback::feedbackWarning(
        inputId = id,
        show = TRUE,
        text = "String could not be coerced to formula."
      )
      req(FALSE)
    }

    a <- try(attributes(terms(formula)))
    if (inherits(a, "try-error")) {
      shinyFeedback::feedbackWarning(inputId = id, show = TRUE, text = "Invalid formula.")
      req(FALSE)
    }

    ok_lhs <-
      a$response > 0L &&
      is.call(lhs <- a$variables[[1L + a$response]]) &&
      lhs[[1L]] == "cbind" &&
      length(lhs) == 3L
    shinyFeedback::feedbackWarning(
      inputId = id,
      show = !ok_lhs,
      text = "Left hand side must be a call to `cbind` with 2 arguments."
    )
    req(ok_lhs)

    ok_rhs <-
      a$intercept == 1L &&
      is.null(a$offset) &&
      length(a$term.labels) < 2L
    shinyFeedback::feedbackWarning(
      inputId = id,
      show = !ok_rhs,
      text = "Right hand side must be 1 or have exactly one term."
    )
    req(ok_rhs)

    attr(formula, "group") <- length(a$term.labels) == 1L
    formula
  }
  make_data <- function(input, id) {
    is_rds <- tools::file_ext(input[[id]]$name) == "rds"
    shinyFeedback::feedbackWarning(
      inputId = id,
      show = !is_rds,
      text = "Invalid file extension."
    )
    req(is_rds)
    data <- readRDS(input[[id]]$datapath)
    is_data_frame <- is.data.frame(data)
    shinyFeedback::feedbackWarning(
      inputId = id,
      show = !is_data_frame,
      text = "Supplied R object must be a data frame."
    )
    req(is_data_frame)
    data
  }
  make_frame <- function(formula, data) {
    lhs <- formula[[2L]]
    rhs <- formula[[3L]]
    call <- call("data.frame",
      g = call("as.factor", rhs),
      cbind1 = lhs[[2L]],
      cbind2 = lhs[[3L]]
    )
    frame <- try(eval(call, data, baseenv()))
    if (inherits(frame, "try-error")) {
      validate(conditionMessage(attr(frame, "condition")))
      req(FALSE)
    }
    attr(frame, "names_original") <- vapply(call[-1L], deparse, "")
    frame
  }
  standardize_missing <- function(x) {
    if (is.double(x)) {
      x[!is.finite(x)] <- NA
    }
    x
  }

  formula <- reactive({
    req(input$formula, input$data)
    make_formula(input, "formula")
  })
  data <- reactive({
    req(input$formula, input$data)
    make_data(input, "data")
  })
  frame <- reactive({
    ff <- make_frame(formula(), data())
    nf <- as.list(attr(ff, "names_original"))
    names(ff) <- names(nf) <- c("ts", "time", "x")

    if (nrow(ff) == 0L) {
      validate(sprintf("Frame constructed from `%s` has zero rows.", deparse(formula())))
      req(FALSE)
    }
    ff$ts <- droplevels(ff$ts)
    if (nlevels(ff$ts) == 0L) {
      validate(sprintf("`%s` must have at least one nonempty level.", nf$ts))
      req(FALSE)
    }
    if (any(table(ff$ts) < 2L)) {
      validate(sprintf("Each level of `%s` must have at least two rows of data.", nf$ts))
      req(FALSE)
    }
    ff <- ff[!is.na(ff$ts), , drop = FALSE]
    if (inherits(ff$time, c("Date", "POSIXt"))) {
      ff$time <- julian(ff$time)
    }
    if (!(is.numeric(ff$time) && all(is.finite(ff$time)))) {
      validate(sprintf("`%s` must be a finite numeric, Date, or POSIXt vector.", nf$time))
      req(FALSE)
    }
    if (nlevels(ff$ts) == 1L) {
      if (!all(diff(ff$time)) > 0) {
        validate(sprintf("`%s` must be increasing.", nf$time))
        req(FALSE)
      }
    } else {
      if (!all(tapply(ff$time, ff$ts, function(x) all(diff(x) > 0)))) {
        validate(sprintf("`%s` must be increasing in each level of `%s`.", nf$time, nf$ts))
        req(FALSE)
      }
    }
    if (!(is.numeric(ff$x) && all(ff$x[!is.na(ff$x)] >= 0))) {
      validate(sprintf("`%s` must be a non-negative numeric vector.", nf$time))
      req(FALSE)
    }
    ff$x <- standardize_missing(ff$x)

    ff <- ff[order(ff$ts), , drop = FALSE]
    row.names(ff) <- NULL
    ff
  })
  K <- reactive(nlevels(frame()$ts))
  n <- reactive(c(tapply(!is.na(frame()$x), frame()$ts, sum)))

  formula_windows <- reactive({
    req(input$formula_windows, input$data_windows)
    make_formula(input, "formula_windows")
  })
  data_windows <- reactive({
    req(input$formula_windows, input$data_windows)
    make_data(input, "data_windows")
  })
  frame_windows <- reactive({
    ff <- make_frame(formula_windows(), data_windows())
    nf <- as.list(attr(ff, "names_original"))
    names(ff) <- names(nf) <- c("ts", "start", "end")

    for (s in c("start", "end")) {
      if (inherits(ff[[s]], c("Date", "POSIXt"))) {
        ff[[s]] <- julian(ff[[s]])
      }
      if (!is.numeric(ff[[s]])) {
        validate(sprintf("`%s` must be a numeric, Date, or POSIXt vector.", nf[[s]]))
        req(FALSE)
      }
      ff[[s]] <- standardize_missing(ff[[s]])
    }

    ff$ts <- factor(ff$ts, levels = levels(frame()$ts))
    ff <- ff[complete.cases(ff), , drop = FALSE]
    ff <- ff[do.call(order, unname(ff)), , drop = FALSE]

    forward <- ff$start < ff$end
    if (!all(forward)) {
      validate(sprintf("Fitting windows (%s, %s] must satisfy `%s < %s`.", nf$start, nf$end, nf$start, nf$end))
      req(FALSE)
    }
    disjoint <- c(by(ff[c("start", "end")], droplevels(ff$ts), function(d) {
      nrow(d) < 2L || all(d$start[-1L] >= d$end[-nrow(d)])
    }))
    if (!all(disjoint)) {
      group <- attr(formula_windows(), "group")
      s <- if (group) sprintf(" in each level of `%s`", nf$ts) else ""
      validate(sprintf("Fitting windows (%s, %s] must be disjoint%s.", nf$start, nf$end, s))
      req(FALSE)
    }

    ff <- ff[order(ff$ts), , drop = FALSE]
    row.names(ff) <- NULL
    ff
  })
  accumulator <- reactiveVal()
  observeEvent(frame_windows(), {
    accumulator(frame_windows())
  })

  output$ui_tabset_ts_select <- renderUI(
    selectInput(
      inputId = "ts",
      label = "Choose a time series.",
      choices = `names<-`(seq_len(K()), levels(frame()$ts))
    )
  )
  observeEvent(frame(), {
    updateTabsetPanel(
      inputId = "tabset_ts",
      selected = if (K() > 1L) "tabset_ts_select" else "tabset_ts_null"
    )
  })
  output$ui_tabset_spar <- renderUI({
    f <- function(i) {
      tabPanelBody(sprintf("tabset_spar_%d", i),
        sliderInput(
          inputId = sprintf("spar_%d", i),
          label = "Cubic spline smoothing parameter:",
          value = 0.66,
          min = 0,
          max = 1,
          step = 0.01
        )
      )
    }
    args <- list(
      id = "tabset_spar",
      type = "hidden",
      tabPanelBody("tabset_spar_null")
    )
    do.call(tabsetPanel, c(args, lapply(seq_len(K()), f)))
  })
  observeEvent(input$ts, {
    updateTabsetPanel(
      inputId = "tabset_spar",
      selected = if (n()[as.integer(input$ts)] >= 4L) sprintf("tabset_spar_%s", input$ts) else "tabset_spar_null"
    )
  })
  output$ui_tabset_plot <- renderUI({
    f <- function(i) {
      tabPanelBody(sprintf("tabset_plot_%d", i),
        uiOutput(sprintf("plot_%d_caption", i)),
        plotOutput(sprintf("plot_%d", i))
      )
    }
    args <- list(
      id = "tabset_plot",
      type = "hidden",
      tabPanelBody("tabset_plot_null")
    )
    do.call(tabsetPanel, c(args, lapply(seq_len(K()), f)))
  })
  observeEvent(input$ts, {
    updateTabsetPanel(
      inputId = "tabset_plot",
      selected = sprintf("tabset_plot_%s", input$ts)
    )
  })
  observeEvent(frame(), {
    f <- function(i) {
      tx <- frame()[unclass(frame()$ts) == i, c("time", "x"), drop = FALSE]
      tx$y <- 1 + tx$x / c(NA, diff(tx$time))
      output[[sprintf("plot_%d_caption", i)]] <- renderUI(HTML("<b>Caption!</b>"))
      output[[sprintf("plot_%d", i)]] <- renderPlot({
        par(
          mar = c(3, 5, 0.5, 1),
          xaxs = "i",
          yaxs = "i",
          las = 1,
          pch = 16,
          cex.lab = 1.2
        )
        plot.new()
        plot.window(xlim = range(tx$time), ylim = range(tx$y, na.rm = TRUE), log = "y")
        points(y ~ time, data = tx, col = "#BBBBBBA8")
        spar <- input[[sprintf("spar_%d", i)]]
        if (n()[i] >= 4L && spar > 0) {
          ss <- smooth.spline(tx$time[!is.na(tx$y)], tx$y[!is.na(tx$y)], spar = spar)
          pp <- predict(ss, seq.int(tx$time[1L], tx$time[nrow(tx)], length.out = 151L))
          lines(y ~ x, data = pp, col = "#004488CC", lwd = 3)
        }
        box()
        epigrowthfit:::Daxis(
          minor = list(mgp = c(3, 0.25, 0), tcl = -0.2, lwd.ticks = 1, gap.axis = 0, cex.axis = 1),
          major = list(mgp = c(3, 1.5, 0), tcl = 0, lwd.ticks = 0, gap.axis = 0, cex.axis = 1.2)
        )
        axis(side = 2, mgp = c(3, 0.7, 0))
        title(ylab = "number of cases per day", line = 4)
      })
      invisible(NULL)
    }
    lapply(seq_len(K()), f)
  })
  output$ui_tabset_points <- renderUI({
    f <- function(i) {
      tabPanelBody(sprintf("tabset_points_%d", i),
        uiOutput(sprintf("points_%d_caption", i)),
        dataTableOutput(sprintf("points_%d", i))
      )
    }
    args <- list(
      id = "tabset_points",
      type = "hidden",
      tabPanelBody("tabset_points_null")
    )
    do.call(tabsetPanel, c(args, lapply(seq_len(K()), f)))
  })
  observeEvent(input$ts, {
    updateTabsetPanel(
      inputId = "tabset_points",
      selected = sprintf("tabset_points_%s", input$ts)
    )
  })
  observeEvent(frame(), {
    tx_split <- split(frame()[c("time", "x")], frame()$ts)
    f <- function(i) {
      output[[sprintf("points_%d_caption", i)]] <- renderUI(HTML("<b>Observations:</b>"))
      output[[sprintf("points_%d", i)]] <- renderDataTable(tx_split[[i]],
        options = list(
          dom = "t",
          paging = FALSE,
          scrollY = "300px"
        )
      )
      invisible(NULL)
    }
    lapply(seq_len(K()), f)
  })
  output$ui_tabset_windows <- renderUI({
    f <- function(i) {
      tabPanelBody(sprintf("tabset_windows_%d", i),
        uiOutput(sprintf("windows_%d_caption", i)),
        dataTableOutput(sprintf("windows_%d", i))
      )
    }
    args <- list(
      id = "tabset_windows",
      type = "hidden",
      tabPanelBody("tabset_windows_null")
    )
    do.call(tabsetPanel, c(args, lapply(seq_len(K()), f)))
  })
  observeEvent(input$ts, {
    updateTabsetPanel(
      inputId = "tabset_windows",
      selected = sprintf("tabset_windows_%s", input$ts)
    )
  })
  observeEvent(accumulator(), {
    se_split <- split(accumulator()[c("start", "end")], accumulator()$ts)
    f <- function(i) {
      output[[sprintf("windows_%d_caption", i)]] <- renderUI(HTML("<b>Fitting windows:</b>"))
      output[[sprintf("windows_%d", i)]] <- renderDataTable(se_split[[i]],
        options = list(
          dom = "t",
          paging = FALSE,
          scrollY = "300px"
        )
      )
    }
    lapply(seq_len(K()), f)
  })
}

shinyApp(ui = ui, server = server)
