## UI for the information about the pipeline
infoUI <- function(datasets) {

  selector <- selectInput("select_pipeline",
                          label="Select data set:",
                          choices=datasets, width="90%")
  if(length(datasets) < 2L) {
    selector <- hidden(selector)
  }

  column(width=12,
    box(title="Overview", width=6,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        column(width=12,
        fluidRow(selector),
        fluidRow(
           tableOutput("project_overview")
    ))),
    box(title="Contrasts", width=6,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        fluidRow(
           tableOutput("contrasts_overview")
    )),
    box(title="Covariates:", width=12,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        column(width=11,
        fluidRow(
             DTOutput("covariates")
    ))),
    box(title="Session info:", width=12,
        solidHeader=TRUE, status="primary", collapsible=TRUE,
        column(width=11,
        fluidRow(
             verbatimTextOutput("session_info")
    )))
        
  )

}


## UI for the help system
helpUI <- function() {
  help_dir <- system.file("seapiper_manual/", package="seaPiper")

  help_files <- list.files(help_dir, full.names = TRUE)

  ret <- lapply(help_files,
                function(f) {
                  title <- gsub("\\.md$", "", gsub("_", " ", basename(f)))
                  box(width=12, collapsible = TRUE, title = title,
                      solidHeader=TRUE, includeMarkdown(f), status="primary")
                })

  do.call(tagList, ret)
}

## add the icon and a dropdown menu to page head
.pipeline_dashboard_header <- function(title) {

  dashboardHeader(title=img(src="icons/piper_horiz.png", alt="[seaPiper]"),
    tags$li(class="dropdown", h4(title, style="font-size:22px;color:white;padding-right:20px;"))
  )

}


.pipeline_dashboard_sidebar <- function(debug_panel=FALSE) {
  menus <- list(
         tipify(menuItem("Gene browser",  tabName = "gene_browser", icon = icon("dna")),
                "Browse genes and view gene expression", placement="right"),
         tipify(menuItem("Volcano plots",  tabName = "volcano_plots", icon = icon("mountain",
                                                                                  class="fa-rotate-180")),
                "Browse genes and view gene expression", placement="right"),
         tipify(menuItem("Tmod browser",  tabName = "tmod_browser", icon = icon("project-diagram")),
                "Browse gene set enrichments", placement="right"),
         tipify(menuItem("Disco plots",   tabName = "disco", icon = icon("chart-line")),
                "Compare contrasts", placement="right"),
         tipify(menuItem("Panel plots",   tabName = "panel_plot", icon = icon("grip-vertical")),
                "Overview of gene set enrichments", placement="right"),
         tipify(menuItem("PCA",           tabName = "pca", icon = icon("cube")),
                "Principal component analysis", placement="right"),
         tipify(menuItem("Workflow Info", tabName = "pip_info", icon = icon("info-circle")),
                "View workflow parameters", placement="right"),
         menuItem("Help",          tabName = "help", icon = icon("question-circle"))
         )

  if(debug_panel) {
    menus <- c(menus, list(menuItem("|DeBuG|", tabName = "os", icon = icon("bug"))))
  }


  dashboardSidebar(
       tags$head(
                 tags$link(rel = "stylesheet", type = "text/css", href = "css/seapiper.css"),
                 tags$link(rel="icon", type="image/png", sizes="32x32", href="icons/favicon-32x32.png"),
                 tags$link(rel="icon", type="image/png", sizes="16x16", href="icons/favicon-16x16.png"),
                 tags$link(rel="shortcut icon", type="image/x-icon", sizes="16x16", href="icons/favicon.ico")),
         # Setting id makes input$tabs give the tabName of currently-selected tab
       do.call(sidebarMenu, c(list(id="navid"), menus))
  )

}


## prepare the actual tabs UI
.pipeline_dashboard_body <- function(data, title, debug_panel=FALSE) {

  npip <- length(data[[1]])

  pipelines <- names(data[[1]])

  cntr_titles <- data[["cntr_titles"]]
  covar       <- data[["covar"]]

      t1 <- tabItem("gene_browser",
         box(title=p("Gene table ", icon("question-circle")),
                     width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserTableUI("geneT", cntr_titles)),
         box(title="Gene info",  width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, geneBrowserPlotUI("geneP", contrasts=TRUE)),
        useShinyjs()
      )
    t10 <- tabItem("volcano_plots", 
         box(title="Volcano plots", width=12, status="primary", 
             collapsible=TRUE,
             solidHeader=TRUE, volcanoUI("volcano", pipelines))
         )
 
    t2 <- tabItem("tmod_browser",
       box(title="Gene set enrichment overview", width=12, status="primary",
           collapsible=TRUE,
           solidHeader=TRUE, tmodBrowserTableUI("tmodT", cntr_titles, upset_pane=TRUE)),
       box(title="Evidence plot", width=12, status="primary",
           collapsible=TRUE,
           solidHeader=TRUE, tmodBrowserPlotUI("tmodP"))
    )
    t3 <- tabItem("disco",
       box(title="Discordance / concordance plots", width=12, status="primary",
       height="800px", solidHeader=TRUE, discoUI("disco", cntr_titles)),
       )
    t4 <- tabItem("panel_plot",
       box(title="Panel plot", width=12, status="primary",
           solidHeader=TRUE, tmodPanelPlotUI("panelP", pipelines)))
    t5 <- tabItem("pca",
       box(title="Principal Component Analysis", width=12, status="primary",
       solidHeader=TRUE, pcaUI("pca", pipelines)),
       useShinyjs()
       )
    t6 <- tabItem("pip_info", 
         infoUI(pipelines))
         
    t7 <- tabItem("help", helpUI())
    t8 <- tabItem("os", fluidPage(h1("test"), tableOutput("debugTab")))


    if(debug_panel) {
      tbit <- tabItems( t1, t10, t2, t3, t4, t5, t6, t7, t8)
    } else {
      tbit <- tabItems( t1, t10, t2, t3, t4, t5, t6, t7)
    }

  dashboardBody(
    tbit,
    style="min-height:2500px;"
  )
}
