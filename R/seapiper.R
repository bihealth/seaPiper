

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
  help_dir <- system.file("seapiper_manual/", package="Rseasnap")

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
#' @param x an object of the class [`seapiper_ds`], for example created with
#' [new_seapiper_dataset()].
#' @param annot_default default annotation which will be used if annotation
#' is not specified in `x`
#' @param tmod_dbs_default default tmod database list, which will be used if the
#' databases are not specified in the `x`
#' @param primary_id name of the column in the annotion data frame
#'        which corresponds to the primary gene identifier (including the
#'        row names of the contrasts results in the cntr object)
#' @param title Name of the pipeline to display
#' @param only_data return the processed data and exit
#' @param lfc_col,pval_col column names for log fold change and p-value
##' @param save_memory if TRUE, then large objects (contrasts, annotations)
##' will be used as disk.frame object. Slower, but more memory efficient.
#' @param debug_panel show a debugging panel
#' @importFrom disk.frame as.disk.frame
#' @importFrom purrr %>%
#' @importFrom shiny renderImage tags img icon imageOutput includeMarkdown addResourcePath
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
#' @importFrom purrr imap map map_chr map_lgl transpose
#' @importFrom tibble rownames_to_column
#' @importFrom DT datatable DTOutput renderDT 
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
seapiper <- function(x, title="Workflow output explorer", 
                             annot_default=NULL, 
                             tmod_dbs_default=NULL,
                             tmod_map_default=NULL,
                             primary_id="PrimaryID",
                             only_data=FALSE, 
                             save_memory=FALSE,
                             lfc_col="log2FoldChange", pval_col="padj",
                             debug_panel=FALSE) {

  ## for debugging purposes
  env <- environment()  # can use globalenv(), parent.frame(), etc
  #env <- .GlobalEnv

  ## x can be a seasnap_ds or a list of seasnap_ds. In this first case, we
  ## change everything into a list.
  if(is(x, "seapiper_ds")) {
    x <- list(default=x)
  }

  ## check that all objects in x are of type seapiper_ds
  stopifnot(all(map_lgl(x, ~ is(.x, "seapiper_ds"))))

  x <- .prepare_data(x, primary_id, annot=annot_default, 
                                       tmod_dbs=tmod_dbs_default,
                                       tmod_map=tmod_map_default,
                                       pval_col=pval_col, lfc_col=lfc_col)
  #data <- .prepare_data(pip, primary_id, annot=annot, cntr=cntr, tmod_res=tmod_res, 
  #                      tmod_dbs=tmod_dbs, save_memory=save_memory)

  if(only_data) { return(x) }

  options(spinner.color="#47336F")
  options(spinner.type=6)

  addResourcePath("icons", system.file("icons", package="Rseasnap"))
  addResourcePath("css",   system.file("css", package="Rseasnap"))

  theme_set(theme_bw())
  thematic_shiny(font="auto")

  ## Prepare the UI
  header  <- .pipeline_dashboard_header(title)     
  sidebar <- .pipeline_dashboard_sidebar(debug_panel=debug_panel)
  body    <- .pipeline_dashboard_body(x, title, debug_panel=debug_panel)
  ui      <- dashboardPage(header, sidebar, body, skin="purple", title=title)

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
                   

    if(is.null(x[["config"]])) {
      output$project_overview   <- renderText({ "No configuration defined" })
    } else {
      output$project_overview   <- renderTable({ project_overview_table(x[["config"]][[ds]], title) })
    }

    if(is.null(x[["config"]])) {
      output$project_overview   <- renderText({ "No configuration defined" })
    } else {
      output$contrasts_overview <- renderTable({ contrasts_overview_table(x[["config"]][[ds]]) }) 
    }

    output$covariates         <- renderDT({ 
      covariate_table(x[["covar"]][[ds]]) 
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

    geneBrowserTableServer("geneT", x[["cntr"]], x[["annot"]], 
      annot_linkout=x[["annot_linkout"]],
      gene_id=gene_id)

    geneBrowserPlotServer("geneP", gene_id, covar=x[["covar"]], 
                          exprs=x[["exprs"]], annot=x[["annot"]], 
                          annot_linkout=x[["annot_linkout"]],
                          cntr=x[["cntr"]],
                          exprs_label = "Regularized log transformed expression (rlog)"
    )

    volcanoServer("volcano", x[["cntr"]], annot=x[["annot"]], gene_id=gene_id)

    discoServer("disco", x[["cntr"]], x[["annot"]], gene_id=gene_id)

    ## whenever gene_id$id changes, navigate to the gene browser
    observeEvent(gene_id$id, {
      updateTabItems(session, "navid", "gene_browser")
    })

    if(debug_panel) {
      message("Using debug panel")
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

 
    tmodBrowserPlotServer("tmodP", gs_id, 
                                      tmod_dbs=x[["tmod_dbs"]], 
                                      cntr    =x[["cntr"]], 
                                      tmod_map=x[["tmod_map"]], 
                                      tmod_gl =x[["tmod_gl"]], 
                                      annot   =x[["annot"]],
                                      tmod_res=x[["tmod_res"]],
                                      gene_id=gene_id)
 
 
    tmodBrowserTableServer("tmodT", x[["tmod_res"]], gs_id=gs_id, 
                           multilevel=TRUE, tmod_dbs=x[["tmod_dbs"]])
    tmodPanelPlotServer("panelP", cntr    =x[["cntr"]], 
                                  tmod_res=x[["tmod_res"]],
                                  tmod_dbs=x[["tmod_dbs"]], 
                                  tmod_map=x[["tmod_map"]], 
                                  annot   =x[["annot"]],
                                  gs_id=gs_id)
 
    if(!is.null(x[["pca"]])) { pcaServer("pca", x[["pca"]], x[["covar"]]) }

    observeEvent(gs_id$id, { updateTabItems(session, "navid", "tmod_browser") })
 
  }

  shinyApp(ui, server, options=list(display.mode="showcase"))
}


## create a table with the overview of the project
project_overview_table <- function(config, title) {
  
  tmp1 <- data.frame(
                 c("Organism", "Taxon ID", "Formula", "Contrasts", "Tmod dbs", 
                   "Low count filter (absolute)", "Low count filter (samples)"),
                 c(config$organism$name, 
                   config$organism$taxon,
                   config$experiment$design_formula,
                   length(config$contrasts$contrast_list),
                   paste(map_chr(config$tmod$databases, `[[`, "title"), collapse=", "),
                   sprintf("Removed genes with < %d total counts", config$filter$low_counts),
                   sprintf("Kept genes with at least %d counts in at least %d samples",
                           config$filter$min_counts, config$filter$min_count_n)
                   ))
  colnames(tmp1) <- c("Title:", title)
  tmp1
}

contrasts_overview_table <- function(config) {

  data.frame(ID=map_chr(config$contrasts$contrast_list, `[[`, "ID"),
             Title=map_chr(config$contrasts$contrast_list, `[[`, "title")
             )


}

## create a datatable with covariates
covariate_table <- function(covar) {

  tmp <- covar[ , !colnames(covar) %in% c("filename", "md5") ]

  datatable(tmp, extensions="Buttons",
              options=list(dom="Bfrtip", scrollX=TRUE, 
                           buttons=c("copy", "csv", "excel")))

}


.prep_annot_linkout <- function(annot, config) {

  ret <- list()
  cn <- colnames(annot)

  if("ENSEMBL" %in% cn) {
    ret$ENSEMBL <- "https://www.ensembl.org/id/%s/"
  }
  if("ENSEMBLID" %in% cn) {
    ret$ENSEMBLID <- "https://www.ensembl.org/id/%s/"
  }
              
  if(config$organism$name == "human" & "SYMBOL" %in% cn) {
    ret$SYMBOL <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
  }

  if("ENTREZ" %in% cn) {
    ret$ENTREZ <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("ENTREZID" %in% cn) {
    ret$ENTREZID <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("REFSEQID" %in% cn) {
    ret$REFSEQID <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  if("REFSEQ" %in% cn) {
    ret$REFSEQ <- "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
  }

  ret
}
