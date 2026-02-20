#' Build seaPiperData from basic inputs
#'
#' Create a `seaPiperData` object from contrasts, annotations, and expression
#' data without requiring a pipeline object.
#' @details
#' If you pass a custom `config`, it must include:
#' `organism` (list with `name`, `taxon`),
#' `experiment` (list with `design_formula`),
#' `contrasts` (list with `contrast_list`, each having `ID` and `title`),
#' `tmod` (list with `databases` and `sort_by`),
#' `filter` (list with `low_counts`, `min_counts`, `min_count_n`),
#' `dataset_title` (string, defaults to "default").
#' @param cntr named list of contrast data frames
#' @param annot annotation data frame
#' @param exprs expression matrix (genes x samples)
#' @param primary_id primary gene identifier column name
#' @param covar optional covariate data frame (samples x variables)
#' @param config optional configuration list; if NULL, defaults are created
#' @param title optional dataset title stored in `config$dataset_title`
#' @export
make_seapiperdata <- function(cntr, annot, exprs,
                              primary_id="PrimaryID",
                              covar=NULL,
                              config=NULL,
                              title=NULL) {
  stopifnot(is.list(cntr))
  stopifnot(is.data.frame(annot))
  stopifnot(is.matrix(exprs) || is.data.frame(exprs))

  if(is.data.frame(exprs)) {
    exprs <- as.matrix(exprs)
  }

  if(is.null(colnames(exprs))) {
    stop("`exprs` must have column names (sample identifiers)")
  }

  if(is.null(rownames(exprs))) {
    stop("`exprs` must have row names (gene identifiers)")
  }

  if(!primary_id %in% colnames(annot)) {
    if(!is.null(rownames(annot))) {
      annot[[primary_id]] <- rownames(annot)
    } else {
      stop(sprintf("`annot` must contain column %s or have row names", primary_id))
    }
  }

  if(is.null(names(cntr))) {
    names(cntr) <- paste0("contrast_", seq_along(cntr))
  }

  cntr <- map(cntr, ~ {
    .x <- as.data.frame(.x)
    if(!primary_id %in% colnames(.x)) {
      if(!is.null(rownames(.x))) {
        .x[[primary_id]] <- rownames(.x)
      } else {
        stop(sprintf("contrast data frames must contain %s or have row names", primary_id))
      }
    }
    .x
  })

  if(is.null(covar)) {
    covar <- data.frame(ID=colnames(exprs), stringsAsFactors=FALSE)
    rownames(covar) <- covar[["ID"]]
  }

  if(is.null(config)) {
    config <- list(
      organism   = list(name="unknown", taxon="unknown"),
      experiment = list(design_formula="not_set"),
      contrasts  = list(contrast_list=map(names(cntr), ~ list(ID=.x, title=.x))),
      tmod       = list(databases=list(), sort_by=character(0)),
      filter     = list(low_counts=NA_integer_, min_counts=NA_integer_, min_count_n=NA_integer_),
      dataset_title = if(is.null(title)) "default" else title
    )
  } else if(is.null(config$dataset_title)) {
    if(!is.null(title)) {
      config$dataset_title <- title
    } else {
      config$dataset_title <- "default"
    }
  }

  cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
  names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")
  cntr_titles <- cntr_titles[ cntr_titles %in% names(cntr) ]

  mtx  <- t(exprs)
  vars <- apply(mtx, 2, var)
  sel  <- vars > 1e-26
  pca  <- prcomp(mtx[ , sel], scale.=TRUE)$x

  data <- list(
    annot= list(default=annot),
    cntr = list(default=cntr),
    tmod_res = list(default=list()),
    tmod_dbs = list(default=list()),
    tmod_map = list(default=NULL),
    tmod_gl  = list(default=list()),
    config   = list(default=config),
    covar    = list(default=covar),
    annot_linkout = list(default=.prep_annot_linkout(annot, config)),
    dbs      = list(default=character(0)),
    sorting  = list(default=character(0)),
    rld      = list(default=exprs),
    pca      = list(default=pca),
    cntr_titles = list(default=cntr_titles)
  )

  structure(data, class="seaPiperData")
}


#' Convert Rseasnap pipeline output into seaPiperData
#'
#' Prepare seaPiper-ready data from a Rseasnap pipeline (or list of pipelines).
#' @param pip pipeline returned by `load_de_pipeline`
#' @param annot annotation returned by `get_annot`
#' @param cntr contrasts list returned by `get_contrasts`
#' @param tmod_dbs tmod databases returned by `get_tmod_dbs`
#' @param tmod_res tmod results returned by `get_tmod_res`
#' @param primary_id name of the column in the annotion data frame
#'        which corresponds to the primary gene identifier (including the
#'        row names of the contrasts results in the cntr object)
#' @importFrom Rseasnap load_de_pipeline 
#' @importFrom Rseasnap get_tmod_res get_tmod_dbs get_tmod_mapping get_config
#' @importFrom Rseasnap get_covariates get_object get_annot get_contrasts 
#' @export
rseasnap_to_seapiperdata <- function(pip, primary_id="PrimaryID",
                                     annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL) {

  ## pip can be a pipeline or a list of pipelines. In this first case, we
  ## change everything into a list.
  if(is(pip, "seasnap_DE_pipeline")) {
    pip <- list(default=pip)

    if(!is.null(annot))    { annot    <- list(default=annot)    }
    if(!is.null(cntr))     { cntr     <- list(default=cntr)     }
    if(!is.null(tmod_res)) { tmod_res <- list(default=tmod_res) }
    if(!is.null(tmod_dbs)) { tmod_dbs <- list(default=tmod_dbs) }
  }

  message("preparing...")
  if(is.null(annot))    { annot <- list() }
  if(is.null(cntr))     { cntr <- list() }
  if(is.null(tmod_res)) { tmod_res <- list() }
  if(is.null(tmod_dbs)) { tmod_dbs <- list() }

  data <- imap(pip, ~ {
    .pip <- .x
    .id  <- .y
    .prepare_data_single_pipeline(.id, .pip, primary_id, annot, cntr, tmod_res, tmod_dbs)
  })

  structure(transpose(data), class="seaPiperData")
}
