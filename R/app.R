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
#' @param pip pipeline returned by `load_de_pipeline`
#' @param annot annotation returned by `get_annot`
#' @param cntr contrasts list returned by `get_contrasts`
#' @param tmod_dbs tmod databases returned by `get_tmod_dbs`
#' @param tmod_res tmod results returned by `get_tmod_res`
#' @param primary_id name of the column in the annotion data frame
#'        which corresponds to the primary gene identifier (including the
#'        row names of the contrasts results in the cntr object)
#' @param title Name of the pipeline to display
#' @param only_data return the processed data and exit
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
#' @importFrom Rseasnap load_de_pipeline 
#' @importFrom Rseasnap get_tmod_res get_tmod_dbs get_tmod_mapping get_config
#' @importFrom Rseasnap get_covariates get_object get_annot get_contrasts 
#' @examples
#' if(interactive()) {
#'   example_dir <- system.file("extdata/example_pipeline", package="Rseasnap")
#'   conf_f      <- file.path(example_dir, "DE_config.yaml")
#'   pip         <- load_de_pipeline(config_file = conf_f)
#'   pipeline_browser(pip)
#' }
#' @export
seapiper <- function(pip, title="Workflow output explorer", 
                             annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL,
                             primary_id="PrimaryID",
                             only_data=FALSE, 
                             debug_panel=FALSE) {
  env <- environment()  # can use globalenv(), parent.frame(), etc
  env <- .GlobalEnv

  ## pip can be a pipeline or a list of pipelines. In this first case, we
  ## change everything into a list.
  if(is(pip, "seasnap_DE_pipeline")) {
    pip <- list(default=pip)

    if(!is.null(annot))    { annot    <- list(default=annot)    }
    if(!is.null(cntr))     { cntr     <- list(default=cntr)     }
    if(!is.null(tmod_res)) { tmod_res <- list(default=tmod_res) }
    if(!is.null(tmod_dbs)) { tmod_dbs <- list(default=tmod_dbs) }
  }

  data <- .prepare_data(pip, primary_id, annot=annot, cntr=cntr, tmod_res=tmod_res, 
                        tmod_dbs=tmod_dbs)
  if(only_data) { return(data) }

  options(spinner.color="#47336F")
  options(spinner.type=6)

  addResourcePath("icons", system.file("icons", package="Rseasnap"))
  addResourcePath("css",   system.file("css", package="Rseasnap"))

  theme_set(theme_bw())
  thematic_shiny(font="auto")

  ## Prepare the UI
  header  <- .pipeline_dashboard_header(title)     
  sidebar <- .pipeline_dashboard_sidebar(debug_panel=debug_panel)
  body    <- .pipeline_dashboard_body(data, title, debug_panel=debug_panel)
  ui <- dashboardPage(header, sidebar, body, skin="purple", title=title)

  #   theme = bs_theme(primary = "#47336F", secondary = "#C6B3EB", 
  #                           font_scale = NULL, 
  #                           `enable-shadows` = TRUE, 
  #                           bootswatch = "united"),


  server <- function(input, output, session) {

    ## pipeline browser specific functions
    observeEvent(input$select_pipeline, {
      ds <- input$select_pipeline
      if(!isTruthy(ds)) {
        return(NULL)
      }
                   
    output$project_overview   <- renderTable({ 
      project_overview_table(data[["config"]][[ds]], title) 
    })
    output$contrasts_overview <- renderTable({ 
      contrasts_overview_table(data[["config"]][[ds]]) })
    output$covariates         <- renderDT({ 
      covariate_table(data[["covar"]][[ds]]) 
      })
    output$session_info <- renderPrint(
                                       sessionInfo()
                                       )
    })

    ## this reactive value holds the id of the selected gene, however the
    ## selection has been done
    gene_id <- reactiveValues()

    ## this is reactive value for gene sets
    gs_id   <- reactiveValues()

    geneBrowserTableServer("geneT", data[["cntr"]], data[["annot"]], 
      annot_linkout=data[["annot_linkout"]],
      gene_id=gene_id)

    tmodBrowserPlotServer("tmodP", gs_id, 
                                      tmod_dbs=data[["tmod_dbs"]], 
                                      cntr    =data[["cntr"]], 
                                      tmod_map=data[["tmod_map"]], 
                                      tmod_gl =data[["tmod_gl"]], 
                                      annot   =data[["annot"]],
                                      tmod_res=data[["tmod_res"]],
                                      gene_id=gene_id)
 
    discoServer("disco", data[["cntr"]], data[["annot"]], gene_id=gene_id)
 
    tmodBrowserTableServer("tmodT", data[["tmod_res"]], gs_id=gs_id, 
                           multilevel=TRUE, tmod_dbs=data[["tmod_dbs"]])
    tmodPanelPlotServer("panelP", cntr    =data[["cntr"]], 
                                  tmod_res=data[["tmod_res"]],
                                  tmod_dbs=data[["tmod_dbs"]], 
                                  tmod_map=data[["tmod_map"]], 
                                  annot   =data[["annot"]],
                                  gs_id=gs_id)
 
    pcaServer("pca", data[["pca"]], data[["covar"]])
    volcanoServer("volcano", data[["cntr"]], annot=data[["annot"]], gene_id=gene_id)
    
    observeEvent(gs_id$id, { 
      updateTabItems(session, "navid", "tmod_browser")
    })
 
    ## combine events selecting a gene from gene browser and from disco
    observeEvent(gene_id$id, {
      updateTabItems(session, "navid", "gene_browser")
    })
 
    geneBrowserPlotServer("geneP", gene_id, covar=data[["covar"]], 
                          exprs=data[["rld"]], annot=data[["annot"]], 
                          annot_linkout=data[["annot_linkout"]],
                          cntr=data[["cntr"]],
                          exprs_label = "Regularized log transformed expression (rlog)"
    )

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
