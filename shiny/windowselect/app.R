library("shiny")
source("utils.R")

ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  titlePanel(
    title = div(HTML("<b>epigrowthfit</b>: Fitting window selection tool")),
    windowTitle = "epigrowthfit: Fitting window selection tool"
  ),
  tabsetPanel(
    id = "download_tabset",
    type = "hidden",
    tabPanelBody(
      value = "download_tab_null"
    ),
    tabPanelBody(
      value = "download_tab_button",
      downloadButton(
        outputId = "download",
        label = HTML("Download session results as <tt>.rds</tt>")
      ),
      br(),
      br()
    )
  ),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "side_tabset",
        type = "tabs",
        tabPanel(
          title = "Upload",
          value = "side_tab_upload",
          br(),
          h4("Time series"),
          h6(HTML("<i>required</i>")),
          textInput(
            inputId = "formula",
            label = "formula",
            placeholder = "e.g., cbind(time, x) ~ ts"
            #value = "cbind(time, x) ~ ts"
          ),
          fileInput(
            inputId = "data",
            label = HTML("data frame [<tt>.rds</tt>]"),
            multiple = FALSE,
            accept = ".rds"
          ),
          h4("Existing fitting windows"),
          h6(HTML("<i>optional</i>")),
          textInput(
            inputId = "formula_windows",
            label = "formula",
            placeholder = "e.g., cbind(start, end) ~ ts"
            #value = "cbind(start, end) ~ ts"
          ),
          fileInput(
            inputId = "data_windows",
            label = HTML("data frame [<tt>.rds</tt>]"),
            multiple = FALSE,
            accept = ".rds"
          )
        ),
        tabPanel(
          title = "Display",
          value = "side_tab_display",
          br(),
          tabsetPanel(
            id = "ts_tabset",
            type = "hidden",
            tabPanelBody(
              value = "ts_tab_null"
            ),
            tabPanelBody(
              value = "ts_tab_select",
              uiOutput(
                outputId = "ui_ts_tab_select"
              )
            )
          ),
          selectInput(
            inputId = "timeas",
            label = "Displayed time format",
            choices = c("Date", "numeric")
          ),
          uiOutput(
            outputId = "ui_xlim_tabset"
          ),
          uiOutput(
            outputId = "ui_logylim_tabset"
          ),
          uiOutput(
            outputId = "ui_spar_tabset"
          )
        )
      )
    ),
    mainPanel(
      uiOutput(
        outputId = "ui_main_tabset"
      )
    )
  )
)

server <- function(input, output, session) {
  make_formula <- function(input, id) {
    formula <- try(as.formula(str2lang(input[[id]]), env = .GlobalEnv), silent = TRUE)
    if (inherits(formula, "try-error")) {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "String could not be coerced to formula.")
      req(FALSE)
    }

    a <- try(attributes(terms(formula)), silent = TRUE)
    if (inherits(a, "try-error")) {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "Invalid formula.")
      req(FALSE)
    }

    ok_lhs <-
      a$response > 0L &&
      is.call(lhs <- a$variables[[1L + a$response]]) &&
      lhs[[1L]] == "cbind" &&
      length(lhs) == 3L
    if (!ok_lhs) {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "Left hand side must be a call to `cbind` with 2 arguments.")
      req(FALSE)
    }

    ok_rhs <-
      a$intercept == 1L &&
      is.null(a$offset) &&
      length(a$term.labels) < 2L
    if (!ok_rhs) {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "Right hand side must be 1 or have exactly one term.")
      req(FALSE)
    }

    attr(formula, "group") <- length(a$term.labels) == 1L
    formula
  }
  make_data <- function(input, id) {
    if (tools::file_ext(input[[id]]$name) != "rds") {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "Invalid file extension.")
      req(FALSE)
    }
    data <- readRDS(input[[id]]$datapath)
    if (!is.data.frame(data)) {
      shinyFeedback::showFeedbackWarning(inputId = id, text = "Supplied R object must be a data frame.")
      req(is_data_frame)
    }
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
    frame <- try(eval(call, data, baseenv()), silent = TRUE)
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
    shinyFeedback::hideFeedback(inputId = "formula")
    req(input$formula)
    make_formula(input, "formula")
  })
  data <- reactive({
    shinyFeedback::hideFeedback(inputId = "data")
    req(input$data)
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
    ff$y <- NA_real_
    split(ff$y, ff$ts) <- c(by(ff[c("time", "x")], ff$ts, function(d) c(NA_real_, 1 + d$x[-1L] / diff(d$time)), simplify = FALSE))

    ff <- ff[order(ff$ts), , drop = FALSE]
    row.names(ff) <- NULL
    attr(ff, "names_original") <- unlist(nf)
    ff
  })
  levels_ts <- reactive(levels(frame()$ts))
  nlevels_ts <- reactive(length(levels_ts()))
  n <- reactive(c(tapply(!is.na(frame()$x), frame()$ts, sum)))
  limits <- reactive({
    f <- function(y, f) if (all(argna <- is.na(y))) 1 else f(y[!argna])
    xmin <- frame()$time[!duplicated(frame()$ts)]
    xmax <- frame()$time[!duplicated(frame()$ts, fromLast = TRUE)]
    ymin <- c(tapply(frame()$y, frame()$ts, function(y) f(y, min)))
    ymax <- c(tapply(frame()$y, frame()$ts, function(y) f(y, max)))
    data.frame(xmin, xmax, ymin, ymax)
  })

  formula_windows <- reactive({
    if (!isTruthy(input$formula_windows)) {
      return(NULL)
    }
    shinyFeedback::hideFeedback(inputId = "formula_windows")
    make_formula(input, "formula_windows")
  })
  data_windows <- reactive({
    if (!isTruthy(input$data_windows)) {
      return(NULL)
    }
    shinyFeedback::hideFeedback(inputId = "data_windows")
    make_data(input, "data_windows")
  })
  frame_windows <- reactive({
    if (is.null(formula_windows()) || is.null(data_windows())) {
      ff <- data.frame(
        ts = factor(levels = levels_ts()),
        start = numeric(),
        end = numeric()
      )
      attr(ff, "names_original") <- `names<-`(names(ff), names(ff))
      return(ff)
    }

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

    ff$ts <- factor(ff$ts, levels = levels_ts())
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
    attr(ff, "names_original") <- unlist(nf)
    ff
  })
  accumulator <- reactiveVal()
  observeEvent(frame_windows(), {
    accumulator(frame_windows())
  })
  N <- reactive(c(table(accumulator()$ts)))

  show_side_tab_display <- reactiveVal(FALSE)
  observeEvent(frame(), show_side_tab_display(TRUE))
  observeEvent(show_side_tab_display(), {
    if (show_side_tab_display()) {
      showTab(
        inputId = "side_tabset",
        target = "side_tab_display",
      )
    } else {
      hideTab(
        inputId = "side_tabset",
        target = "side_tab_display",
      )
    }
  })
  output$ui_ts_tab_select <- renderUI(
    selectInput(
      inputId = "ts",
      label = "Displayed time series",
      choices = `names<-`(seq_len(nlevels_ts()), levels_ts())
    )
  )
  outputOptions(output, "ui_ts_tab_select", suspendWhenHidden = FALSE)
  observeEvent(frame(), {
    updateTabsetPanel(
      inputId = "ts_tabset",
      selected = if (nlevels_ts() > 1L) "ts_tab_select" else "ts_tab_null"
    )
  })
  output$ui_xlim_tabset <- renderUI({
    args <- list(
      id = "xlim_tabset",
      type = "hidden",
      tabPanelBody(
        value = "xlim_tab_null"
      )
    )
    f <- function(i) {
      tabPanelBody(
        value = sprintf("xlim_tab_%d_Date", i),
        dateRangeInput(
          inputId = sprintf("xlim_%d_Date", i),
          label = "Horizontal axis limits",
          start = .Date(-1),
          end = .Date(1),
          min = .Date(-1),
          max = .Date(1)
        )
      )
    }
    g <- function(i) {
      tabPanelBody(
        value = sprintf("xlim_tab_%d_numeric", i),
        sliderInput(
          inputId = sprintf("xlim_%d_numeric", i),
          label = "Horizontal axis limits:",
          value = c(-1, 1),
          min = -1,
          max = 1
        )
      )
    }
    do.call(tabsetPanel, c(args, lapply(seq_len(nlevels_ts()), f), lapply(seq_len(nlevels_ts()), g)))
  })
  outputOptions(output, "ui_xlim_tabset", suspendWhenHidden = FALSE)
  output$ui_logylim_tabset <- renderUI({
    args <- list(
      id = "logylim_tabset",
      type = "hidden",
      tabPanelBody(
        value = "logylim_tab_null"
      )
    )
    f <- function(i) {
      tabPanelBody(
        value = sprintf("logylim_tab_%d", i),
        sliderInput(
          inputId = sprintf("logylim_%d", i),
          label = HTML("Vertical axis limits, log<sub>10</sub> scale"),
          value = c(-1, 1),
          min = -1,
          max = 1,
          step = 0.1
        )
      )
    }
    do.call(tabsetPanel, c(args, lapply(seq_len(nlevels_ts()), f)))
  })
  outputOptions(output, "ui_logylim_tabset", suspendWhenHidden = FALSE)
  output$ui_spar_tabset <- renderUI({
    args <- list(
      id = "spar_tabset",
      type = "hidden",
      tabPanelBody(
        value = "spar_tab_null"
      )
    )
    f <- function(i) {
      tabPanelBody(
        value = sprintf("spar_tab_%d", i),
        sliderInput(
          inputId = sprintf("spar_%d", i),
          label = "Smoothing parameter",
          value = 0.66,
          min = 0,
          max = 1,
          step = 0.01
        )
      )
    }
    do.call(tabsetPanel, c(args, lapply(seq_len(nlevels_ts()), f)))
  })
  outputOptions(output, "ui_spar_tabset", suspendWhenHidden = FALSE)
  observeEvent(list(input$timeas, input$ts), {
    updateTabsetPanel(
      inputId = "xlim_tabset",
      selected = sprintf("xlim_tab_%s_%s", input$ts, input$timeas)
    )
    updateTabsetPanel(
      inputId = "logylim_tabset",
      selected = sprintf("logylim_tab_%s", input$ts)
    )
    updateTabsetPanel(
      inputId = "spar_tabset",
      selected = if (n()[[as.integer(input$ts)]] >= 4L) sprintf("spar_tab_%s", input$ts) else "spar_tab_null"
    )
  })
  observeEvent(limits(), {
    f <- function(i) {
      xval <- c(limits()$xmin[i], limits()$xmax[i])
      xlim <- xval + c(-1, 1) * 0.1 * diff(xval)
      updateDateRangeInput(
        inputId = sprintf("xlim_%d_Date", i),
        start = .Date(xval[1L]),
        end = .Date(xval[2L]),
        min = .Date(xlim[1L]),
        max = .Date(xlim[2L])
      )
      updateSliderInput(
        inputId = sprintf("xlim_%d_numeric", i),
        value = xval,
        min = xlim[1L],
        max = xlim[2L]
      )
      logyval <- log10(c(limits()$ymin[[i]], limits()$ymax[[i]]))
      logylim <- c(min(-1, floor(logyval[1L])), ceiling(logyval[2L]))
      updateSliderInput(
        inputId = sprintf("logylim_%d", i),
        value = logyval,
        min = logylim[1L],
        max = logylim[2L]
      )
    }
    lapply(seq_len(nlevels_ts()), f)
  })
  observeEvent(lapply(sprintf("xlim_%d_Date", seq_len(nlevels_ts())), function(id) input[[id]]), {
    req(input$ts)
    i <- as.integer(input$ts)
    updateSliderInput(
      inputId = sprintf("xlim_%d_numeric", i),
      value = julian(input[[sprintf("xlim_%d_Date", i)]])
    )
  })
  observeEvent(lapply(sprintf("xlim_%d_numeric", seq_len(nlevels_ts())), function(id) input[[id]]), {
    req(input$ts)
    i <- as.integer(input$ts)
    updateSliderInput(
      inputId = sprintf("xlim_%d_Date", i),
      value = .Date(input[[sprintf("xlim_%d_numeric", i)]])
    )
  })
  output$ui_main_tabset <- renderUI({
    args <- list(
      id = "main_tabset",
      type = "hidden",
      tabPanelBody(
        value = "main_tab_null"
      )
    )
    f <- function(i) {
      tabPanelBody(
        value = sprintf("main_tab_%d", i),
        fluidRow(
          uiOutput(
            outputId = sprintf("plot_%d_caption", i)
          ),
          plotOutput(
            outputId = sprintf("plot_%d", i),
            # click = clickOpts(
            #   id = sprintf("click_%d", i),
            #   clip = TRUE
            # ),
            dblclick = dblclickOpts(
              id = sprintf("dblclick_%d", i),
              clip = TRUE,
              delay = 400
            ),
            # hover = hoverOpts(
            #   id = sprintf("hover_%d", i),
            #   clip = TRUE,
            #   delay = 300,
            #   delayType = "debounce",
            #   nullOutside = TRUE
            # ),
            brush = brushOpts(
              id = sprintf("brush_%d", i),
              fill = "auto",
              stroke = "auto",
              opacity = 0.25,
              delay = 1000,
              delayType = "debounce",
              clip = TRUE,
              direction = "xy",
              resetOnNew = TRUE
            )
          )
        ),
        br(),
        br(),
        fluidRow(
          column(6,
            uiOutput(outputId = sprintf("points_%d_caption", i)),
            dataTableOutput(outputId = sprintf("points_%d", i))
          ),
          column(6,
            uiOutput(outputId = sprintf("windows_%d_caption", i)),
            dataTableOutput(outputId = sprintf("windows_%d", i))
          )
        )
      )
    }
    do.call(tabsetPanel, c(args, lapply(seq_len(nlevels_ts()), f)))
  })
  outputOptions(output, "ui_main_tabset", suspendWhenHidden = FALSE, priority = 100)
  observeEvent(input$ts, {
    updateTabsetPanel(
      inputId = "main_tabset",
      selected = sprintf("main_tab_%s", input$ts)
    )
  })
  observeEvent(list(input$timeas, frame(), accumulator()), {
    req(input$timeas, frame(), accumulator())
    index1 <- split(seq_len(nrow(frame())), frame()$ts)
    index2 <- split(seq_len(nrow(accumulator())), accumulator()$ts)
    f <- function(i) {
      dp1 <- frame()[index1[[i]], c("time", "y"), drop = FALSE]
      dp2 <- accumulator()[index2[[i]], c("start", "end"), drop = FALSE]

      output[[sprintf("plot_%d_caption", i)]] <- renderUI(HTML(paste0(
        "<small>",
        "To select a fitting window, click and drag the pointer over the plot region. ",
        "<br>",
        "To deselect a fitting window, double click anywhere inside of it.",
        "</small>"
      )))
      output[[sprintf("plot_%d", i)]] <- renderPlot({
        par(
          mar = c(3.5, 5, 0.5, 1),
          xaxs = "i",
          yaxs = "i",
          las = 1,
          pch = 16,
          cex.lab = 1.2
        )
        xlim <- input[[sprintf("xlim_%d_numeric", i)]]
        if (is.na(xlim[1L])) {
          if (is.na(xlim[2L])) {
            xlim <- rep_len(0, 2L)
          } else {
            xlim <- rep_len(xlim[2L], 2L)
          }
        }
        if (is.na(xlim[2L])) {
          xlim[2L] <- xlim[1L]
        }
        ylim <- 10^(input[[sprintf("logylim_%d", i)]])
        plot.new()
        plot.window(xlim = xlim, ylim = ylim, log = "y")
        usr <- par("usr")
        if (nrow(dp2) > 0L) {
          rect(
            xleft = dp2$start,
            xright = dp2$end,
            ybottom = 10^usr[3L],
            ytop = 10^usr[4L],
            col = "#0044881A",
            border = NA
          )
        }
        points(y ~ time, data = dp1, col = "#BBBBBBA8")
        spar <- input[[sprintf("spar_%d", i)]]
        argna <- is.na(dp1$y)
        if (sum(!argna) >= 4L && spar > 0) {
          ss <- smooth.spline(dp1$time[!argna], dp1$y[!argna], spar = spar)
          xx <- seq.int(
            from = max(usr[1L], limits()$xmin[i]),
            to = min(usr[2L], limits()$xmax[i]),
            length.out = 151L
          )
          lines(y ~ x, data = predict(ss, xx), col = "#004488CC", lwd = 3)
        }
        box()
        if (input$timeas == "Date") {
          Daxis(
            minor = list(mgp = c(3, 0.25, 0), tcl = -0.2, lwd.ticks = 1, gap.axis = 0, cex.axis = 1),
            major = list(mgp = c(3, 1.5, 0), tcl = 0, lwd.ticks = 0, gap.axis = 0, cex.axis = 1.2)
          )
        } else {
          axis(side = 1, mgp = c(3, 0.7, 0))
          title(xlab = "time", line = 2.5)
        }
        axis(side = 2, mgp = c(3, 0.7, 0))
        title(ylab = "1 + (number of cases per day)", line = 4)
      })

      data_table_options <- list(dom = "t", paging = FALSE, scrollY = "200px")

      dt1 <- frame()[index1[[i]], c("time", "x"), drop = FALSE]
      if (input$timeas == "Date") {
        dt1$time <- .Date(dt1$time)
      }
      output[[sprintf("points_%d_caption", i)]] <- renderUI(HTML("<b>Observations:</b>"))
      output[[sprintf("points_%d", i)]] <- renderDataTable(dt1, options = data_table_options)

      dt2 <- accumulator()[index2[[i]], c("start", "end"), drop = FALSE]
      if (input$timeas == "Date") {
        dt2[] <- lapply(dt2, .Date)
      }
      output[[sprintf("windows_%d_caption", i)]] <- renderUI(HTML("<b>Fitting windows:</b>"))
      output[[sprintf("windows_%d", i)]] <- renderDataTable(dt2, options = data_table_options)
      invisible(NULL)
    }
    lapply(seq_len(nlevels_ts()), f)
  })

  observeEvent(lapply(sprintf("brush_%d", seq_len(nlevels_ts())), function(id) input[[id]]), {
    req(input$ts)
    i <- as.integer(input$ts)
    brush <- input[[sprintf("brush_%d", i)]]
    req(brush, brush$xmin < brush$xmax)
    time <- frame()$time[unclass(frame()$ts) == i]
    req(brush$xmin < time[length(time)], brush$xmax >= time[2L])
    if (brush$xmin < time[2L]) {
      start <- time[1L]
    } else {
      start <- time[which.max(time > brush$xmin) - 1L]
    }
    if (brush$xmax >= time[length(time)]) {
      end <- time[length(time)]
    } else {
      end <- time[which.max(time > brush$xmax) - 1L]
    }
    req(start < end)
    newrow <- data.frame(ts = factor(levels_ts()[i], levels = levels_ts()), start, end)
    res <- rbind(accumulator(), newrow)
    res <- res[do.call(order, unname(res)), , drop = FALSE]
    if (N()[[i]] > 0L) {
      k <- which(unclass(res$ts) == i)
      req(all(res$start[k[-1L]] >= res$end[k[-length(k)]]))
    }
    accumulator(res)
  })

  observeEvent(lapply(sprintf("dblclick_%d", seq_len(nlevels_ts())), function(id) input[[id]]), {
    req(input$ts)
    i <- as.integer(input$ts)
    usr <- input[[sprintf("dblclick_%d", i)]]
    req(usr)
    deselected <-
      unclass(accumulator()$ts) == i &
      accumulator()$start < usr$x &
      accumulator()$end >= usr$x
    if (any(deselected)) {
      accumulator(accumulator()[!deselected, , drop = FALSE])
    }
  })

  observeEvent(frame(), {
    updateTabsetPanel(
      inputId = "download_tabset",
      selected = "download_tab_button"
    )
  })
  output$download <- downloadHandler(
    filename = function() {
      old <- tools::file_path_sans_ext(input$data$name)
      sprintf("%s_shinyout.rds", old)
    },
    content = function(file) {
      res <- list(
        formula = cbind(time, x) ~ ts,
        data = frame()[c("ts", "time", "x")],
        formula_windows = cbind(start, end) ~ ts,
        data_windows = accumulator()[c("ts", "start", "end")]
      )
      environment(res$formula) <- environment(res$formula_windows) <- .GlobalEnv
      if (input$timeas == "Date") {
        res$data$time <- .Date(res$data$time)
        res$data_windows$start <- .Date(res$data_windows$start)
        res$data_windows$end <- .Date(res$data_windows$start)
      }
      saveRDS(res, file = file)
    }
  )
}

shinyApp(ui = ui, server = server)
