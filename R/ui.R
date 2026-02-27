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

  dashboardHeader(title=img(src="icons/piper_horiz.svg", alt="[seaPiper]", height="40px"),
    tags$li(class="dropdown", h4(title, style="font-size:22px;color:white;padding-right:20px;"))
  )

}


.pipeline_dashboard_sidebar <- function(features, debug_panel=FALSE) {
  menus <- list()

  if(isTRUE(features$gene_browser)) {
    menus <- c(menus, list(
         tipify(menuItem("Gene browser",  tabName = "gene_browser", icon = icon("dna")),
                "Browse genes and view gene expression", placement="right")))
  }

  if(isTRUE(features$volcano)) {
    menus <- c(menus, list(
         tipify(menuItem("Volcano plots",  tabName = "volcano_plots", icon = icon("mountain",
                                                                                  class="fa-rotate-180")),
                "Browse genes and view gene expression", placement="right")))
  }

  if(isTRUE(features$heatmap)) {
    menus <- c(menus, list(
         tipify(menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
                "View expression heatmaps for selected genes", placement="right")))
  }

  if(isTRUE(features$disco)) {
    menus <- c(menus, list(
         tipify(menuItem("Disco plots",   tabName = "disco", icon = icon("chart-line")),
                "Compare contrasts", placement="right")))
  }

  if(isTRUE(features$tmod)) {
    menus <- c(menus, list(
         tipify(menuItem("Tmod browser",  tabName = "tmod_browser", icon = icon("project-diagram")),
                "Browse gene set enrichments", placement="right")))
  }

  if(isTRUE(features$tmod_panel)) {
    menus <- c(menus, list(
         tipify(menuItem("Panel plots",   tabName = "panel_plot", icon = icon("grip-vertical")),
                "Overview of gene set enrichments", placement="right")))
  }

  if(isTRUE(features$pca)) {
    menus <- c(menus, list(
         tipify(menuItem("PCA",           tabName = "pca", icon = icon("cube")),
                "Principal component analysis", placement="right")))
  }

  if(isTRUE(features$covar)) {
    menus <- c(menus, list(
         tipify(menuItem("Color palettes", tabName = "color_palettes", icon = icon("palette")),
                "Configure color palettes", placement="right")))
  }

  if(isTRUE(features$file_export)) {
    menus <- c(menus, list(
         tipify(menuItem("Export data",   tabName = "file_export", icon = icon("download")),
                "Export displayed data objects", placement="right")))
  }

  if(isTRUE(features$info)) {
    menus <- c(menus, list(
         tipify(menuItem("Workflow Info", tabName = "pip_info", icon = icon("info-circle")),
                "View workflow parameters", placement="right")))
  }

  menus <- c(menus, list(menuItem("Help", tabName = "help", icon = icon("question-circle"))))

  if(debug_panel) {
    menus <- c(menus, list(menuItem("|DeBuG|", tabName = "os", icon = icon("bug"))))
  }


  dashboardSidebar(
       tags$head(
                 tags$link(rel = "stylesheet", type = "text/css", href = "css/seapiper.css"),
                 tags$link(rel="icon", type="image/svg+xml", href="icons/favicon.svg"),
                 tags$link(rel="icon", type="image/png", sizes="32x32", href="icons/favicon-32x32.png"),
                 tags$link(rel="icon", type="image/png", sizes="16x16", href="icons/favicon-16x16.png"),
                 tags$link(rel="shortcut icon", type="image/x-icon", sizes="16x16", href="icons/favicon.ico")),
         # Setting id makes input$tabs give the tabName of currently-selected tab
       do.call(sidebarMenu, c(list(id="navid"), menus))
  )

}


## prepare the actual tabs UI
.pipeline_dashboard_body <- function(data, title, features, export_objects=NULL, debug_panel=FALSE) {

  pipelines <- .get_dataset_ids(data)
  if(length(pipelines) == 0) {
    pipelines <- "default"
  }

  pipeline_titles <- pipelines
  if(is.list(data[["config"]]) && length(data[["config"]]) > 0) {
    pipeline_titles <- vapply(pipelines, function(ds) {
      cfg <- data[["config"]][[ds]]
      if(is.null(cfg) || is.null(cfg$dataset_title)) {
        return(ds)
      }
      as.character(cfg$dataset_title)
    }, character(1))
  }
  pipeline_choices <- pipelines
  names(pipeline_choices) <- pipeline_titles

  cntr_titles <- data[["cntr_titles"]]
  if(is.null(cntr_titles) && is.list(data[["cntr"]]) && length(data[["cntr"]]) > 0) {
    cntr_titles <- imap(data[["cntr"]], ~ {
      if(is.null(.x)) {
        return(character(0))
      }
      ids <- names(.x)
      if(is.null(ids)) {
        ids <- paste0("contrast_", seq_along(.x))
      }
      names(ids) <- ids
      ids
    })
    if(length(cntr_titles) == 1L) {
      cntr_titles <- cntr_titles[[1]]
    }
  }

    tabs <- list()

    if(isTRUE(features$gene_browser)) {
      tabs <- c(tabs, list(
        tabItem("gene_browser",
          box(title=p("Gene table ", icon("question-circle")),
                      width=12, status="primary", 
              collapsible=TRUE,
              solidHeader=TRUE, geneBrowserTableUI(id="geneT", cntr_titles=cntr_titles)),
          box(title="Gene info",  width=12, status="primary", 
              collapsible=TRUE,
              solidHeader=TRUE, geneBrowserPlotUI(id="geneP", contrasts=TRUE)),
         useShinyjs()
        )
      ))
    }

    if(isTRUE(features$volcano)) {
      tabs <- c(tabs, list(
        tabItem("volcano_plots", 
          box(title="Volcano plots", width=12, status="primary", 
              collapsible=TRUE,
              solidHeader=TRUE, volcanoUI("volcano", pipeline_choices))
        )
      ))
    }

    if(isTRUE(features$heatmap)) {
      tabs <- c(tabs, list(
        tabItem("heatmap",
          box(title="Expression heatmap", width=12, status="primary",
              collapsible=TRUE,
              solidHeader=TRUE, bioshmods::heatmapUI("heatmap"))
        )
      ))
    }

    if(isTRUE(features$disco)) {
      tabs <- c(tabs, list(
        tabItem("disco",
          box(title="Discordance / concordance plots", width=12, status="primary",
          height="800px", solidHeader=TRUE, discoUI("disco", cntr_titles))
        )
      ))
    }

    if(isTRUE(features$tmod)) {
      tabs <- c(tabs, list(
        tabItem("tmod_browser",
          box(title="Gene set enrichment overview", width=12, status="primary",
              collapsible=TRUE,
              solidHeader=TRUE, tmodBrowserTableUI("tmodT", cntr_titles, upset_pane=TRUE)),
          box(title="Evidence plot", width=12, status="primary",
              collapsible=TRUE,
              solidHeader=TRUE, tmodBrowserPlotUI("tmodP"))
        )
      ))
    }

    if(isTRUE(features$tmod_panel)) {
      tabs <- c(tabs, list(
        tabItem("panel_plot",
          box(title="Panel plot", width=12, status="primary",
              solidHeader=TRUE, tmodPanelPlotUI("panelP", pipeline_choices)))
      ))
    }

    if(isTRUE(features$pca)) {
      tabs <- c(tabs, list(
        tabItem("pca",
          box(title="Principal Component Analysis", width=12, status="primary",
          solidHeader=TRUE, pcaUI("pca", pipeline_choices)),
          useShinyjs()
        )
      ))
    }

    if(isTRUE(features$covar)) {
      tabs <- c(tabs, list(
        tabItem("color_palettes",
          box(title="Color palettes", width=12, status="primary",
              solidHeader=TRUE, bioshmods::colorPalettesUI("color_palettes"))
        )
      ))
    }

    if(isTRUE(features$file_export) && !is.null(export_objects) && length(export_objects) > 0L) {
      tabs <- c(tabs, list(
        tabItem("file_export",
          box(title="Export data", width=12, status="primary",
          solidHeader=TRUE, fileExportUI("fexp", export_objects))
        )
      ))
    }

    if(isTRUE(features$info)) {
      tabs <- c(tabs, list(
        tabItem("pip_info", infoUI(pipeline_choices))
      ))
    }
         
    tabs <- c(tabs, list(tabItem("help", helpUI())))

    if(debug_panel) {
      tabs <- c(tabs, list(tabItem("os", fluidPage(h1("test"), tableOutput("debugTab")))))
    }

    tbit <- do.call(tabItems, tabs)

  dashboardBody(
    tbit,
    style="min-height:2500px;"
  )
}

## determine dataset IDs from the first available list-like section
.get_dataset_ids <- function(data) {
  candidates <- c("covar", "config", "cntr", "annot", "rld", "pca",
                  "tmod_res", "tmod_dbs", "tmod_map", "tmod_gl", "cntr_titles")
  for(field in candidates) {
    value <- data[[field]]
    if(is.list(value) && length(value) > 0 && !is.null(names(value))) {
      ids <- names(value)
      if(length(ids) > 0) {
        return(ids)
      }
    }
  }
  character(0)
}
