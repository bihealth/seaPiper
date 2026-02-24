#' seaPiper Workflow output explorer
#'
#' seaPiper Workflow output explorer for the sea-snap pipeline
#' 
#' Launch a shiny app for viewing / exploring the results of a differential
#' expression analysis performed in sea-snap.
#'
#' Launching is much faster if the large objects (contrasts, tmod
#' databases) do not need to be loaded each time the pipeline browser is
#' started.
#' @param data a `seaPiperData` object created by `seapiperdata_from_rseasnap`,
#'        `seapiperdata_from_objects`, or `seapiperdata_from_yaml`
#' @param title Name of the pipeline to display
#' @param debug_panel show a debugging panel
#' @importFrom purrr %>%
#' @importFrom shiny renderImage tags img icon imageOutput includeMarkdown
#' @importFrom shiny addResourcePath
#' @importFrom shinydashboard dashboardPage dashboardBody dashboardSidebar dashboardHeader
#' @importFrom shinydashboard box sidebarMenu tabItem tabItems updateTabItems menuItem
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom shiny isTruthy shinyApp
#' @importFrom shiny tableOutput renderTable renderPrint verbatimTextOutput
#' @importFrom shiny selectInput 
#' @importFrom shiny observeEvent reactiveValues
#' @importFrom shiny column fluidPage fluidRow h1 h2 h3 h4 p tagList
#' @importFrom shinyjs useShinyjs hidden
#' @importFrom purrr imap map map_chr transpose
#' @importFrom tibble rownames_to_column
#' @importFrom DT datatable
#' @importFrom DT DTOutput renderDT 
#' @importFrom thematic thematic_shiny
#' @importFrom ggplot2 theme_bw theme_set
#' @importFrom shinyBS tipify
#' @importFrom bioshmods geneBrowserTableServer tmodBrowserPlotServer discoServer 
#' @importFrom bioshmods tmodBrowserTableServer tmodPanelPlotServer pcaServer volcanoServer geneBrowserPlotServer
#' @importFrom bioshmods geneBrowserTableUI geneBrowserPlotUI volcanoUI tmodBrowserTableUI
#' @importFrom bioshmods tmodBrowserPlotUI discoUI tmodPanelPlotUI pcaUI
#' @importFrom utils sessionInfo object.size
#' @examples
#' if(interactive()) {
#'   example_dir <- system.file("extdata/example_pipeline", package="Rseasnap")
#'   conf_f      <- file.path(example_dir, "DE_config.yaml")
#'   pip         <- load_de_pipeline(config_file = conf_f)
#'   data        <- seapiperdata_from_rseasnap(pip)
#'   seapiper(data)
#' }
#' @export
seapiper <- function(data, title="Workflow output explorer", 
                             debug_panel=FALSE) {
  env <- environment()  # can use globalenv(), parent.frame(), etc
  env <- .GlobalEnv

  required_bioshmods <- "0.0.0.9004"
  installed_bioshmods <- as.character(utils::packageVersion("bioshmods"))
  if(utils::compareVersion(installed_bioshmods, required_bioshmods) < 0) {
    stop(
      sprintf(
        paste(
          "seaPiper requires bioshmods >= %s, but %s is installed.",
          "Install/update bioshmods from /home/january/Projects/R/bioshmods/bioshmods."
        ),
        required_bioshmods,
        installed_bioshmods
      )
    )
  }

  if(!inherits(data, "seaPiperData")) {
    stop("`data` must be a seaPiperData object")
  }

  options(spinner.color="#47336F")
  options(spinner.type=6)

  addResourcePath("icons", system.file("icons", package="seaPiper"))
  addResourcePath("css",   system.file("css", package="seaPiper"))

  theme_set(theme_bw())
  thematic_shiny(font="auto")

  ## Prepare the UI
  validation <- validate_seapiperdata(data)
  features <- validation$features

  if(length(validation$missing) > 0) {
    details <- vapply(names(validation$missing), function(name) {
      sprintf("%s (missing: %s)", name, paste(validation$missing[[name]], collapse=", "))
    }, character(1))
    message("seaPiper: disabled modules: ", paste(details, collapse="; "))
  }
  header  <- .pipeline_dashboard_header(title)     
  sidebar <- .pipeline_dashboard_sidebar(features=features, debug_panel=debug_panel)
  body    <- .pipeline_dashboard_body(data, title, features=features, debug_panel=debug_panel)
  ui <- dashboardPage(header, sidebar, body, skin="purple", title=title)

  #   theme = bs_theme(primary = "#47336F", secondary = "#C6B3EB", 
  #                           font_scale = NULL, 
  #                           `enable-shadows` = TRUE, 
  #                           bootswatch = "united"),


  server <- function(input, output, session) {

    ## this reactive value holds the id of the selected gene, however the
    ## selection has been done
    gene_id <- reactiveValues()

    ## this is reactive value for gene sets
    gs_id   <- reactiveValues()

    if(isTRUE(features$info)) {
      .seapiper_server_info(input, output, session, data, title)
    }

    if(isTRUE(features$gene_browser)) {
      .seapiper_server_gene(input, output, session, data, gene_id)
    }

    if(isTRUE(features$tmod)) {
      .seapiper_server_tmod(input, output, session, data, gs_id, gene_id,
                            enable_panel=isTRUE(features$tmod_panel))
    }

    .seapiper_server_misc(input, output, session, data, gene_id,
                          enable_disco=isTRUE(features$disco),
                          enable_volcano=isTRUE(features$volcano),
                          enable_pca=isTRUE(features$pca))

    if(debug_panel) {
      output$debugTab <- renderTable({
        objects <- ls(env)
        sizes <- unlist(lapply(objects, function(x) {
                          object.size(get(x, envir=env, inherits=FALSE))}))
        sizes.mb <- unlist(lapply(objects, function(x) {
                          format(object.size(get(x, envir=env, inherits=FALSE)), units="Mb")}))

        ret <- data.frame(
                          objects=objects,
                          sizes  =sizes,
                          sizes.mb=sizes.mb)
        ret <- ret[ order(-ret$sizes), ]
        ret
      })
    }
  }

  shinyApp(ui, server)
}
