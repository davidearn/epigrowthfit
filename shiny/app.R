library("shiny")

ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  withMathJax(),
  titlePanel(
    title = div(HTML("<b>epigrowthfit</b>: Top level nonlinear models")),
    windowTitle = "epigrowthfit: Top level nonlinear models"
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "curve",
        label = "Nonlinear model of expected cumulative incidence:",
        choices = c(
          exponential = "exponential",
          subexponential = "subexponential",
          Gompertz = "gompertz",
          logistic = "logistic",
          Richards = "richards"
        )
      ),
      uiOutput("mathjax_curve"),
      tabsetPanel(
        id = "par_curve",
        type = "hidden",
        tabPanel("exponential",
          sliderInput(
            inputId = "exponential_log_r",
            label = "\\(\\log(r)\\)",
            value = -3,
            min = -9,
            max = 1,
            step = 0.1
          ),
          uiOutput("exponential_r"),
          sliderInput(
            inputId = "exponential_log_c0",
            label = "\\(\\log(c(0))\\)",
            value = 0,
            min = -5,
            max = 5,
            step = 0.1
          ),
          uiOutput("exponential_c0")
        ),
        tabPanel("subexponential",
          sliderInput(
            inputId = "subexponential_log_alpha",
            label = "\\(\\log(\\alpha)\\)",
            value = -2.5,
            min = -9,
            max = 1,
            step = 0.1
          ),
          uiOutput("subexponential_alpha"),
          sliderInput(
            inputId = "subexponential_log_c0",
            label = "\\(\\log(c(0))\\)",
            value = 0,
            min = -5,
            max = 5,
            step = 0.1
          ),
          uiOutput("subexponential_c0"),
          sliderInput(
            inputId = "subexponential_logit_p",
            label = "\\(\\mathrm{logit}(p)\\)",
            value = 3,
            min = -5,
            max = 5,
            step = 0.1
          ),
          uiOutput("subexponential_p")
        ),
        tabPanel("gompertz",
          sliderInput(
            inputId = "gompertz_log_alpha",
            label = "\\(\\log(\\alpha)\\)",
            value = -2.5,
            min = -9,
            max = 1,
            step = 0.1
          ),
          uiOutput("gompertz_alpha"),
          sliderInput(
            inputId = "gompertz_log_tinfl",
            label = "\\(\\log(t_{\\text{infl}})\\)",
            value = 4.5,
            min = -1,
            max = 9,
            step = 0.1
          ),
          uiOutput("gompertz_tinfl"),
          sliderInput(
            inputId = "gompertz_log_K",
            label = "\\(\\log(K)\\)",
            value = 7,
            min = 0,
            max = 10,
            step = 0.1
          ),
          uiOutput("gompertz_K")
        ),
        tabPanel("logistic",
          sliderInput(
            inputId = "logistic_log_r",
            label = "\\(\\log(r)\\)",
            value = -3,
            min = -9,
            max = 1,
            step = 0.1
          ),
          uiOutput("logistic_r"),
          sliderInput(
            inputId = "logistic_log_tinfl",
            label = "\\(\\log(t_{\\text{infl}})\\)",
            value = 4.5,
            min = -1,
            max = 9,
            step = 0.1
          ),
          uiOutput("logistic_tinfl"),
          sliderInput(
            inputId = "logistic_log_K",
            label = "\\(\\log(K)\\)",
            value = 7,
            min = 0,
            max = 10,
            step = 0.1
          ),
          uiOutput("logistic_K"),
        ),
        tabPanel("richards",
          sliderInput(
            inputId = "richards_log_r",
            label = "\\(\\log(r)\\)",
            value = -3,
            min = -9,
            max = 1,
            step = 0.1
          ),
          uiOutput("richards_r"),
          sliderInput(
            inputId = "richards_log_tinfl",
            label = "\\(\\log(t_{\\text{infl}})\\)",
            value = 4.5,
            min = -1,
            max = 9,
            step = 0.1
          ),
          uiOutput("richards_tinfl"),
          sliderInput(
            inputId = "richards_log_K",
            label = "\\(\\log(K)\\)",
            value = 7,
            min = 0,
            max = 10,
            step = 0.1
          ),
          uiOutput("richards_K"),
          sliderInput(
            inputId = "richards_log_a",
            label = "\\(\\log(a)\\)",
            value = 0,
            min = -5,
            max = 5,
            step = 0.1
          ),
          uiOutput("richards_a"),
        )
      ),
      selectInput(
        inputId = "family",
        label = "Family of observation distributions:",
        choices = c(
          Poisson = "pois",
          `negative binomial` = "nbinom"
        )
      ),
      uiOutput("mathjax_family"),
      tabsetPanel(
        id = "par_family",
        type = "hidden",
        tabPanel("pois"),
        tabPanel("nbinom",
          sliderInput(
            inputId = "nbinom_log_disp",
            label = "\\(\\log(k)\\)",
            value = 4,
            min = -9,
            max = 9,
            step = 0.1
          ),
          uiOutput("nbinom_disp")
        )
      ),
      sliderInput(
        inputId = "range_time",
        label = "Displayed time interval:",
        value = c(0L, 150L),
        min = -300L,
        max = 300L,
        step = 10L
      )
    ),
    mainPanel(
      plotOutput("plot", height = "500px"),
      uiOutput("caption")
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$curve, {
    updateTabsetPanel(inputId = "par_curve", selected = input$curve)
  })
  observeEvent(input$family, {
    updateTabsetPanel(inputId = "par_family", selected = input$family)
  })

  output$mathjax_curve <- renderUI(withMathJax(HTML(
    sprintf("\\[c(t) = %s\\]",
      switch(input$curve,
        exponential = "c(0) \\exp(r t)",
        subexponential = "c(0) [1 + \\alpha (1 - p) c(0)^{-(1 - p)} t]^{1 / (1 - p)}",
        gompertz = "K \\exp(-\\exp(-\\alpha (t - t_{\\text{infl}})))",
        logistic = "\\frac{K}{1 + \\exp(-r (t - t_{\\text{infl}}))}",
        richards = "\\frac{K}{[1 + a \\exp(-a r (t - t_{\\text{infl}}))]^{1 / a}}"
      )
    )
  )))
  output$mathjax_family <- renderUI(withMathJax(HTML(
    sprintf("\\[X(s, t) \\sim %s\\]",
      switch(input$family,
        pois = "\\mathrm{Poisson}(\\lambda(s,t)) \\\\ \\lambda(s,t) = c(t) - c(s)",
        nbinom = "\\mathrm{NegativeBinomial}(\\mu(s,t), k) \\\\ \\mu(s,t) = c(t) - c(s)"
      )
    )
  )))

  output$exponential_r        <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{r = %.6f}\\\\\\)",       exp(input$exponential_log_r)))))
  output$exponential_c0       <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{c(0) = %.6f}\\\\\\)",    exp(input$exponential_log_c0)))))
  output$subexponential_alpha <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{\\alpha = %.6f}\\\\\\)", exp(input$subexponential_log_alpha)))))
  output$subexponential_c0    <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{c(0) = %.6f}\\\\\\)",    exp(input$subexponential_log_c0)))))
  output$subexponential_p     <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{p = %0.6f}\\\\\\)",       plogis(input$subexponential_logit_p)))))
  output$gompertz_alpha       <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{\\alpha = %.6f}\\\\\\)", exp(input$gompertz_log_alpha)))))
  output$gompertz_tinfl       <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{t_{\\text{infl}} = %.6f}\\\\\\)", exp(input$gompertz_log_tinfl)))))
  output$gompertz_K           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{K = %.6f}\\\\\\)",       exp(input$gompertz_log_K)))))
  output$logistic_r           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{r = %.6f}\\\\\\)",       exp(input$logistic_log_r)))))
  output$logistic_tinfl       <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{t_{\\text{infl}} = %.6f}\\\\\\)", exp(input$logistic_log_tinfl)))))
  output$logistic_K           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{K = %.6f}\\\\\\)",       exp(input$logistic_log_K)))))
  output$richards_r           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{r = %.6f}\\\\\\)",       exp(input$richards_log_r)))))
  output$richards_tinfl       <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{t_{\\text{infl}} = %.6f}\\\\\\)", exp(input$richards_log_tinfl)))))
  output$richards_K           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{K = %.6f}\\\\\\)",       exp(input$richards_log_K)))))
  output$richards_a           <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{a = %.6f}\\\\\\)",       exp(input$richards_log_a)))))
  output$nbinom_disp          <- renderUI(withMathJax(HTML(sprintf("\\(\\scriptsize{k = %.6f}\\\\\\)",       exp(input$nbinom_log_disp)))))

  eval_log_curve <- function(time, input) {
    switch(input$curve,
      exponential = {
        input$exponential_log_c0 + exp(input$exponential_log_r) * time
      },
      subexponential = {
        one_minus_p <- plogis(-input$subexponential_logit_p)
        input$subexponential_log_c0 + log1p(one_minus_p * exp(input$subexponential_log_alpha - one_minus_p * input$subexponential_log_c0) * time) / one_minus_p
      },
      gompertz = {
        input$gompertz_log_K - exp(-exp(input$gompertz_log_alpha) * (time - exp(input$gompertz_log_tinfl)))
      },
      logistic = {
        input$logistic_log_K - log1p(exp(-exp(input$logistic_log_r) * (time - exp(input$logistic_log_tinfl))))
      },
      richards = {
        a <- exp(input$richards_log_a)
        input$richards_log_K - log1p(a * exp(-a * exp(input$richards_log_r) * (time - exp(input$richards_log_tinfl)))) / a
      }
    )
  }
  eval_log_r <- function(log_curve, input) {
    switch(input$curve,
      exponential = rep_len(input$exponential_log_r, length(log_curve)),
      subexponential = input$subexponential_log_alpha - plogis(-input$subexponential_logit_p) * log_curve,
      gompertz = input$gompertz_log_alpha + log(input$gompertz_log_K - log_curve),
      logistic = input$logistic_log_r + log1p(-exp(log_curve - input$logistic_log_K)),
      richards = input$richards_log_r + log1p(-exp(exp(input$richards_log_a) * (log_curve - input$richards_log_K)))
    )
  }

  n <- reactive({
    n <- diff(input$range_time)
    shinyFeedback::feedbackWarning(
      inputId = "range_time",
      show = (n == 0L),
      text = "Must have nonzero width."
    )
    req(n > 0L)
    n
  })
  time <- reactive({
    short <- seq.int(input$range_time[1L], input$range_time[2L])
    long <- seq.int(input$range_time[1L], input$range_time[2L], length.out = 151L)
    list(short = short, long = long, long_minus_one = long - 1)
  })

  ## log(c(t))
  log_curve <- reactive(lapply(time(), eval_log_curve, input = input))
  ## c(t)
  curve <- reactive(lapply(log_curve(), exp))
  ## c(t) - c(t - 1)
  diff_curve <- reactive(curve()$long - curve()$long_minus_one)
  ## c'(t) / c(t)
  r <- reactive(exp(eval_log_r(log_curve = log_curve()$long, input = input)))

  x <- reactive(switch(input$family,
    pois = c(NA, rpois(n(), lambda = diff(curve()$short))),
    nbinom = c(NA, rnbinom(n(), mu = diff(curve()$short), size = exp(input$nbinom_log_disp)))
  ))
  c0_plus_cumsum_x <- reactive(c(NA, curve()$long[1L] + cumsum(as.numeric(x())[-1L])))
  zz <- reactive({
    if (n() >= 4L) {
      tmp <- log(c0_plus_cumsum_x()[-1L])
      c(NA, NA, (tmp[-(1:2)] - tmp[-(length(tmp) - 1:0)]) / 2, NA)
    } else {
      rep_len(NA_real_, n() + 1L)
    }
  })

  output$plot <- renderPlot(res = 96, expr = {
    par(
      mfrow = c(3L, 2L),
      mar = c(2, 4, 0.5, 0.5),
      oma = c(2, 0, 2, 0),
      xaxs = "i",
      yaxs = "i",
      mgp = c(3, 0.7, 0),
      las = 1,
      pch = 16,
      cex = 0.9
    )
    col_points <- "#BBBBBBCC"
    col_lines <- c("#BB5566CC", "#DDAA33CC", "#004488CC")
    lwd_lines <- 3

    do_plot <- function(y_points, y_lines, y_log, col_points, col_lines, lwd_lines) {
      if (y_log) {
        log <- "y"
        ymin <- max(exp(-12), min(y_points[y_points > 0],   y_lines[y_lines > 0],   na.rm = TRUE))
        ymax <- min(exp(+12), max(y_points[y_points < Inf], y_lines[y_lines < Inf], na.rm = TRUE))
        ylim <- c(ymin, ymax) * (ymax / ymin)^(c(-1, 1) * 0.04)
      } else {
        log <- ""
        ylim <- c(0, max(y_points[y_points < Inf], y_lines[y_lines < Inf], na.rm = TRUE) * 1.04)
      }
      plot.new()
      plot.window(xlim = input$range_time, ylim = ylim, log = log)
      points(time()$short, y_points, col = col_points)
      lines(time()$long, y_lines, col = col_lines, lwd = lwd_lines)
      box()
      axis(side = 1)
      axis(side = 2)
      invisible(NULL)
    }

    do_plot(y_points = c0_plus_cumsum_x(), y_lines = curve()$long, y_log = FALSE,
            col_points = col_points, col_lines = col_lines[1L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c)(italic(t))))
    title(main = "linear scale", line = 1, xpd = NA)
    do_plot(y_points = c0_plus_cumsum_x(), y_lines = curve()$long, y_log = TRUE,
            col_points = col_points, col_lines = col_lines[1L], lwd_lines = lwd_lines)
    title(main = "logarithmic scale", line = 1, xpd = NA)

    do_plot(y_points = x(), y_lines = diff_curve(), y_log = FALSE,
            col_points = col_points, col_lines = col_lines[2L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c)(italic(t)) - italic(c)(italic(t) - 1)))
    do_plot(y_points = x(), y_lines = diff_curve(), y_log = TRUE,
            col_points = col_points, col_lines = col_lines[2L], lwd_lines = lwd_lines)

    do_plot(y_points = zz(), y_lines = r(), y_log = FALSE,
            col_points = col_points, col_lines = col_lines[3L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c) * "'" * (italic(t)) * " " / " " * italic(c)(italic(t))))
    title(sub = "time", line = 2, xpd = NA)
    do_plot(y_points = zz(), y_lines = r(), y_log = TRUE,
            col_points = col_points, col_lines = col_lines[3L], lwd_lines = lwd_lines)
    title(sub = "time", line = 2, xpd = NA)
  })

  output$caption <- renderUI(withMathJax(HTML(paste0(
    "<b>Row 1:</b> ",
    "Line \\(c(t)\\) shows expected cumulative incidence. ",
    "Points \\((t_{i}, y_{i})\\) show observed cumulative incidence. ",
    "<br/><br/>",
    "<b>Row 2:</b> ",
    "Line \\(c(t) - c(t - 1)\\) shows expected interval incidence. ",
    "Points \\((t_{i}, x_{i})\\) show observed interval incidence. ",
    "<br/><br/>",
    "<b>Row 3:</b> ",
    "Line \\(c'(t) / c(t)\\) shows the predicted per capita growth rate. ",
    "Points \\((t_{i}, r_{i})\\) show a central difference approximation.",
    "<br/><br/>",
    "<b>Notation:</b> \\[\\begin{aligned} ",
    "t_{i} &= ", input$range_time[1L], " + i\\,, \\\\ ",
    "x_{i} &= [\\text{realization of } X(t_{i} - 1,t_{i})]\\,, \\\\ ",
    "y_{i} &= c(t_{0}) + \\sum_{j = 1}^{i} x_{j}\\,, \\\\ ",
    "r_{i} &= \\frac{\\log(y_{i+1}) - \\log(y_{i-1})}{2}\\,. ",
    "\\end{aligned}\\]"
  ))))
}

shinyApp(ui, server)
