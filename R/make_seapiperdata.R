#' Build seaPiperData from basic inputs
#'
#' Create a `seaPiperData` object from contrasts, annotations, and expression
#' data without requiring a pipeline object.
#' @param cntr named list of contrast data frames
#' @param annot annotation data frame
#' @param exprs expression matrix (genes x samples)
#' @param primary_id primary gene identifier column name
#' @param covar optional covariate data frame (samples x variables)
#' @export
make_seapiperdata <- function(cntr, annot, exprs,
                              primary_id="PrimaryID",
                              covar=NULL) {
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

  config <- list(
    organism   = list(name="unknown", taxon="unknown"),
    experiment = list(design_formula="not_set"),
    contrasts  = list(contrast_list=map(names(cntr), ~ list(ID=.x, title=.x))),
    tmod       = list(databases=list(), sort_by=character(0)),
    filter     = list(low_counts=NA_integer_, min_counts=NA_integer_, min_count_n=NA_integer_)
  )

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
