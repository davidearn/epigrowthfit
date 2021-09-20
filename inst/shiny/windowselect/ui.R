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
