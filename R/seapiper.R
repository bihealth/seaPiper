

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


## load missing objects, prepare the data etc.
.prepare_data <- function(pip, primary_id, annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL, save_memory=FALSE) {

  message("preparing...")
  if(is.null(annot))    { annot <- list() }
  if(is.null(cntr))     { cntr <- list() }
  if(is.null(tmod_res)) { tmod_res <- list() }
  if(is.null(tmod_dbs)) { tmod_dbs <- list() }

  data <- imap(pip, ~ {
    .pip <- .x
    .id  <- .y
    .prepare_data_single_pipeline(.id, .pip, primary_id, annot, cntr, tmod_res, tmod_dbs, save_memory=save_memory)
  })

  return(transpose(data))
}

## prepares the data structure for a single pipeline
.prepare_data_single_pipeline <- function(.id, .pip, primary_id, annot, cntr, tmod_res, tmod_dbs,
                                          contrast_cols=c(primary_id, "log2FoldChange", "pvalue", "padj"),
                                          save_memory=FALSE) {
  ret <- list()

  if(is.null(annot[[.id]])) {
    message(sprintf(" * Loading annotation for %s (consider using the annot option to speed this up)", .id))
    ret[["annot"]] <- get_annot(.pip)
  } else {
    ret[["annot"]] <- annot[[.id]]
  }

  if(save_memory) {
    ret[["annot"]] <- as.disk.frame(ret[["annot"]])
  }

  if(is.null(cntr[[.id]])) {
    message(sprintf(" * Loading contrasts for %s (consider using the cntr option to speed this up)", .id))
    ret[["cntr"]] <- get_contrasts(.pip)
  } else {
    ret[["cntr"]] <- cntr[[.id]]
  }

  ret[["cntr"]] <- map(ret[["cntr"]], ~ {
                    ret <- .x %>% rownames_to_column(primary_id) 
                    ret[ , colnames(ret) %in% contrast_cols ]
                                          })
  if(save_memory) {
    ret[["cntr"]] <- map(ret[["cntr"]], ~ as.disk.frame(.x))
  }

  if(is.null(tmod_res[[.id]])) {
    message(sprintf(" * Loading tmod results for %s (consider using the tmod_res option to speed this up)", .id))
    ret[["tmod_res"]] <- get_tmod_res(.pip)
  } else {
    ret[["tmod_res"]] <- tmod_res[[.id]]
  }

  if(is.null(tmod_dbs[[.id]])) {
    message(sprintf(" * Loading tmod databases for %s (consider using the tmod_dbs option to speed this up)", .id))
    ret[["tmod_dbs"]] <- get_tmod_dbs(.pip)
  } else {
    ret[["tmod_dbs"]] <- tmod_dbs[[.id]]
  }

  ## get rid of unnecessary data
  for(i in 1:length(ret[["tmod_dbs"]])) {
    ret[["tmod_dbs"]][[i]][["dbobj"]][["GENES2MODULES"]] <- NULL
  }

  ## we only want the tmod object
  ret[["tmod_dbs"]] <- map(ret[["tmod_dbs"]], ~ .x$dbobj)

  ret[["tmod_map"]] <- get_tmod_mapping(.pip)
  ret[["tmod_gl"]]  <- get_object(.pip, step="tmod", extension="gl.rds", as_list=TRUE)

  ## we only need the order of the genes from one tmod db
  ret[["tmod_gl"]] <- map(ret[["tmod_gl"]], ~  # one for each of contrast
                             map(.[[1]], ~  # one for each sorting type
                                 match(names(.), ret[["annot"]][[primary_id]])))

  ret[["config"]]   <- get_config(.pip)
  ret[["covar"]]    <- get_covariates(.pip)

  ret[["annot_linkout"]] <- .prep_annot_linkout(ret[["annot"]], ret[["config"]])

  ret[["dbs"]]     <- names(tmod_dbs)
  ret[["sorting"]] <- ret[["config"]]$tmod$sort_by

  ret[["rld"]]     <- get_object(.pip, step="DESeq2", extension="rld.blind.rds")
  ret[["rld"]]     <- assay(ret[["rld"]])

  ## prepare the PCA
  mtx <- t(ret[["rld"]])
  vars <- order(apply(mtx, 2, var), decreasing=TRUE)
  sel  <- vars > 1e-26
  ret[["pca"]] <- prcomp(mtx[ , sel], scale.=TRUE)$x
  
  ret[["cntr_titles"]]        <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "ID")
  names(ret[["cntr_titles"]]) <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "title")
  ret[["cntr_titles"]]        <- ret[["cntr_titles"]][ ret[["cntr_titles"]] %in% names(ret[["cntr"]]) ]

  ret
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
#' @param save_memory if TRUE, then large objects (contrasts, annotations)
#' will be used as disk.frame object. Slower, but more memory efficient.
#' @param debug_panel show a debugging panel
#' @importFrom disk.frame as.disk.frame
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
                             save_memory=FALSE,
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
                        tmod_dbs=tmod_dbs, save_memory=save_memory)
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
