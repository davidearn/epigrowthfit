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
        label = div(HTML("Upload an <tt>.rds</tt> file.")),
        multiple = FALSE,
        accept = ".rds"
      ),
      textInput(
        inputId = "formula",
        label = div(HTML("Provide a formula indicating the data structure.")),
        value = "",
        placeholder = "x ~ time | ts"
      )
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  formula <- reactive({
    req(input$formula)
    formula <- try(as.formula(str2lang(input$formula)))
    is_formula <- inherits(formula, "formula")
    shinyFeedback::feedbackWarning(
      inputId = "formula",
      show = !is_formula,
      text = "String could not be coerced to formula."
    )
    req(is_formula)
    a <- attributes(terms(formula))
    ok <-
      a$response > 0L &&
      a$intercept == 1L &&
      is.null(a$offset) &&
      length(a$term.labels) == 1L
    shinyFeedback::feedbackWarning(
      inputId = "formula",
      show = !ok,
      text = "Invalid formula."
    )
    req(ok)
    formula
  })

  data <- reactive({
    req(input$data)
    is_rds <- tools::file_ext(input$data$name) == "rds"
    shinyFeedback::feedbackWarning(
      inputId = "rds",
      show = !is_rds,
      text = "Invalid file extension."
    )
    req(is_rds)
    data <- readRDS(input$data$datapath)
    is_data_frame <- is.data.frame(data)
    shinyFeedback::feedbackWarning(
      inputId = "rds",
      show = !is_data_frame,
      text = "Supplied R object must be a data frame."
    )
    req(is_data_frame)
    data
  })

  frame <- reactive({
    lhs <- formula()[[2L]]
    rhs <- formula()[[3L]]
    l <- list(quote(data.frame), time = lhs, x = rhs)
    is_bar <- function(x) is.call(x) && x[[1L]] == as.name("|")
    if (is_bar <- is.call(rhs) && rhs[[1L]] == as.name("|")) {
      l[c("time", "ts")] <- rhs[2:3]
    }
    frame <- eval(as.call(l), envir = data())
    if (!is_bar) {
      frame$ts <- rep_len(factor(1), nrow(frame))
    }
    nf <- lapply(l[-1L], deparse)

    need(nrow(frame) > 0L, "Supplied data frame has 0 rows.")

    if (inherits(frame$time, "Date")) {
      frame$time <- julian(frame$time, origin = min(frame$time))
    }
    need(is.numeric(frame$time) && all(is.finite(frame$time)),
         sprintf("`%s` must be a finite numeric or Date vector.", nf$time))
    need(is.numeric(frame$x) && frame$x[!is.na(frame$x)] >= 0,
         sprintf("`%s` must be a non-negative numeric vector.", nf$x))


  })




}
