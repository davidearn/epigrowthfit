ui <- fluidPage(
  withMathJax(),
  titlePanel(
    title = div(HTML("<b>epigrowthfit</b>: Top level nonlinear models")),
    windowTitle = "epigrowthfit: Top level nonlinear models"
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "curve",
        label = "Nonlinear model of expected cumulative disease incidence:",
        choices = c(
          exponential = "exponential",
          subexponential = "subexponential",
          Gompertz = "gompertz",
          logistic = "logistic",
          Richards = "richards"
        ),
        selected = "logistic"
      ),
      uiOutput("mathjax_curve"),
      tabsetPanel(
        id = "par_curve",
        type = "hidden",
        tabPanel(
          title = "exponential",
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
        tabPanel(
          title = "subexponential",
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
        tabPanel(
          title = "gompertz",
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
        tabPanel(
          title = "logistic",
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
          uiOutput("logistic_K")
        ),
        tabPanel(
          title = "richards",
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
          uiOutput("richards_a")
        )
      ),
      selectInput(
        inputId = "family",
        label = "Family of observation distributions:",
        choices = c(
          Poisson = "pois",
          `negative binomial` = "nbinom"
        ),
        selected = "nbinom"
      ),
      uiOutput("mathjax_family"),
      tabsetPanel(
        id = "par_family",
        type = "hidden",
        tabPanel(
          title = "pois"
        ),
        tabPanel(
          title = "nbinom",
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
