server <- function(input, output, session) {
  observeEvent(input$curve, {
    updateTabsetPanel(inputId = "par_curve", selected = input$curve)
  })
  observeEvent(input$family, {
    updateTabsetPanel(inputId = "par_family", selected = input$family)
  })
  observeEvent(input$range_time, {
    if (input$range_time[1L] == input$range_time[2L]) {
      updateSliderInput(
        inputId = "range_time",
        value = c(min(290L, input$range_time[1L]), min(300L, input$range_time[1L] + 10L))
      )
    }
  })

  output$mathjax_curve <- renderUI(withMathJax(HTML(
    sprintf("\\[c(t) = %s\\]",
      switch(input$curve,
        exponential = "c(0) \\exp(r t)",
        subexponential = "c(0) \\big(1 + \\tfrac{1 - p}{c(0)^{1 - p}} \\alpha t\\big)^{1 / (1 - p)}",
        gompertz = "\\frac{K}{\\exp(\\exp(-\\alpha (t - t_{\\text{infl}})))}",
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

  n <- reactive(diff(input$range_time) + 1L)
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
                       pois = rpois(n(), lambda = c(curve()$short[1L], diff(curve()$short))),
                       nbinom = rnbinom(n(), mu = c(curve()$short[1L], diff(curve()$short)), size = exp(input$nbinom_log_disp))
  ))
  cumsum_x <- reactive(cumsum(as.numeric(x())))
  zz <- reactive({
    if (n() >= 3L) {
      tmp <- log(cumsum_x())
      c(NA, (tmp[-(1:2)] - tmp[-(n() - 1:0)]) / 2, NA)
    } else {
      rep_len(NA_real_, n())
    }
  })

  output$plot <- renderPlot(res = 96, expr = {
    par(
      mfrow = c(3L, 2L),
      mar = c(2.5, 4.2, 0.1, 0.7),
      oma = c(1, 0, 2, 0),
      xaxs = "i",
      yaxs = "i",
      mgp = c(3, 0.7, 0),
      las = 1,
      pch = 16,
      cex = 0.9
    )
    col_points <- "#BBBBBBA8"
    col_lines <- c("#BB5566CC", "#DDAA33CC", "#004488CC")
    lwd_lines <- 3

    do_plot <- function(y_points, y_lines, y_log, col_points, col_lines, lwd_lines) {
      if (y_log) {
        if (all(y_points == 0, y_lines == 0, na.rm = TRUE)) {
          plot.new()
          return(invisible(NULL))
        }
        log <- "y"
        ymin <- min(y_points[y_points > 0], y_lines[y_lines > 0], na.rm = TRUE)
        ymax <- max(y_points[y_points < Inf], y_lines[y_lines < Inf], na.rm = TRUE)
        if (sum(y_lines > exp(-12), na.rm = TRUE) / sum(y_lines > 0, na.rm = TRUE) > 0.5) {
          ymin <- max(exp(-12), ymin)
        }
        if (sum(y_lines < exp(12), na.rm = TRUE) / sum(y_lines < Inf, na.rm = TRUE) > 0.5) {
          ymax <- min(exp(12), ymax)
        }
        ylim <- exp(log(c(ymin, ymax)) + c(-1, 1) * 0.08 * (log(ymax) - log(ymin)))
        ylim[ylim == 0 | ylim == Inf] <- c(ymin, ymax)[ylim == 0 | ylim == Inf]
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

    do_plot(y_points = cumsum_x(), y_lines = curve()$long, y_log = FALSE,
            col_points = col_points, col_lines = col_lines[1L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c)(italic(t))))
    title(main = "linear", line = 1, xpd = NA)
    do_plot(y_points = cumsum_x(), y_lines = curve()$long, y_log = TRUE,
            col_points = col_points, col_lines = col_lines[1L], lwd_lines = lwd_lines)
    title(main = "logarithmic", line = 1, xpd = NA)

    do_plot(y_points = c(NA, x()[-1L]), y_lines = diff_curve(), y_log = FALSE,
            col_points = col_points, col_lines = col_lines[2L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c)(italic(t)) - italic(c)(italic(t) - 1)))
    do_plot(y_points = c(NA, x()[-1L]), y_lines = diff_curve(), y_log = TRUE,
            col_points = col_points, col_lines = col_lines[2L], lwd_lines = lwd_lines)

    do_plot(y_points = zz(), y_lines = r(), y_log = FALSE,
            col_points = col_points, col_lines = col_lines[3L], lwd_lines = lwd_lines)
    title(ylab = expression(italic(c) * "'" * (italic(t)) * " " / " " * italic(c)(italic(t))))
    title(sub = expression("time, " * italic(t)), line = 2, xpd = NA)
    do_plot(y_points = zz(), y_lines = r(), y_log = TRUE,
            col_points = col_points, col_lines = col_lines[3L], lwd_lines = lwd_lines)
    title(sub = expression("time, " * italic(t)), line = 2, xpd = NA)
  })

  output$caption <- renderUI(withMathJax(HTML(paste0(
    "<b>Row 1:</b> ",
    "Line \\(c(t)\\) shows expected cumulative disease incidence. ",
    "Points \\((t_{i}, y_{i})\\) show simulated cumulative disease incidence. ",
    "<br/><br/>",
    "<b>Row 2:</b> ",
    "Line \\(c(t) - c(t - 1)\\) shows expected disease incidence. ",
    "Points \\((t_{i}, x_{i})\\) show simulated disease incidence. ",
    "<br/><br/>",
    "<b>Row 3:</b> ",
    "Line \\(c'(t) / c(t)\\) shows the predicted per capita growth rate. ",
    "Points \\((t_{i}, r_{i})\\) show a central difference approximation.",
    "<br/><br/>",
    "<b>Notation:</b> \\[\\begin{aligned} ",
    "t_{i} &=", if (input$range_time[1L] == 0L) "" else paste(input$range_time[1L], "+"), " i\\,, \\\\ ",
    "x_{i} &= \\text{realization of} \\begin{cases} ",
    "X(-\\infty,t_{0})\\,, & i = 0\\,, \\\\ ",
    "X(t_{i} - 1,t_{i})\\,, & i = 1,2,\\ldots, ",
    "\\end{cases} \\\\ ",
    "y_{i} &= \\sum_{j = 0}^{i} x_{j}\\,, \\\\ ",
    "r_{i} &= \\frac{\\log(y_{i+1}) - \\log(y_{i-1})}{2}\\,. ",
    "\\end{aligned}\\]"
  ))))
}
