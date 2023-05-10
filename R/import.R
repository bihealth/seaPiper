# check that data frame contains primary ID. If not
# create primary ID column and populate it with
# row names. 
.check_df_primary_id <- function(x, primary_id) {
  if(!is.data.frame(x)) { x <- as.data.frame(x) }
  if(primary_id %in% colnames(x)) { return(x) }

  x[[primary_id]] <- row.names(x)
  x
}

#' Construct a new seapiper data set from individual objects
#'
#' Construct a new seapiper data set from individual objects
#'
#' @param cntr Contrasts: a list of contrasts. If the list is unnamed, then
#' names will be generated (ID1, ID2, ...). Each element of the list must
#' be a data frame or an object that can be converted to a data frame. If
#' the column with primary identifiers (i.e., column with the name
#' "PrimaryID" or whatever is given by the primary_id parameter) is missing, then row names will be
#' saved to that column. If both the column with primary IDs and row names
#' are missing, an error will be raised.
#' @param tmod_res Optional: Results of the gene set enrichment analysis. If the list is unnamed, then
#' names will be generated (ID1, ID2, ...). The sets of names of `tmod_res`
#' and `cntr` must be equal (that is, there should be an element for each
#' contrast). Each element of the list must
#' itself be either a data frame with gene set enrichment results or a named list of data frames
#' with gene set enrichment results
#' @param annot Annotation data frame. If defined, the column indicated by
#' the `primary_id` parameter must be present; alternatively, the rows
#' define the primary gene IDs. If both `primary_id` and row names are
#' missing, an error is raised.
#' @param tmod_dbs Optional: A named list of gene set databases. See "Details".
#' @param tmod_map Optional: an object defining mapping between the gene
#' set databases and the `primary_id` column from the contrasts. See
#' "Details".
#' @param primary_id Name of the column in contrast data frames and in the
#' annotation data frame which holds the primary gene identifier.
#' @param exprs data frame or matrix containing the normalized expression
#' values
#' @return An object of the class seapiper_ds
#' @rdname seapiper_ds
#' @export
# XXX how should we deal with tmod gene sets results? biglist or list of
# dataframes?
new_seapiper_dataset <- function(
                                 cntr,
                                 covar,
                                 exprs      = NULL,
                                 cntr_titles= NULL,
                                 tmod_res   = NULL,
                                 annot      = NULL,
                                 primary_id = "PrimaryID",
                                 tmod_dbs   = NULL,
                                 tmod_map   = NULL

  
  #primary_id, annot, cntr, tmod_res, tmod_dbs,
                                          #contrast_cols=c(primary_id, "log2FoldChange", "pvalue", "padj"),
                                          #save_memory=FALSE
                                 ) {
  ret <- list()

  stopifnot(is.list(cntr))
  stopifnot(is.list(tmod_res))

  if(is.null(names(cntr))) {
    names(cntr) <- paste0("ID", seq_along(cntr))
  }

  message("checking contrasts")
  ## XXX this is unnecessarily slow because we copy the whole object even if it
  ## is OK... I think
  cntr <- imap(cntr, ~ {
    # XXX check for .x being a DF
    .x <- .check_df_primary_id(.x, primary_id)
    .x
  })

  ret$cntr  <- cntr

  if(is.null(cntr_titles)) { cntr_titles <- names(cntr) }

  ret$cntr_titles <- cntr_titles
  names(ret$cntr_titles) <- names(cntr)

  # store the gene expression object
  if(is.matrix(exprs)) {
    message("Converting exprs matrix to data frame")
    exprs <- as.data.frame(exprs)
  }
  ret$exprs       <- exprs

  stopifnot(is.data.frame(covar))
  ret$covar <- covar

  message("checking tmod res")
  if(!is.null(tmod_res)) {
    if(is.null(names(tmod_res))) {
      names(tmod_res) <- paste0("ID", seq_along(tmod_res))
    }

    stopifnot(setequal(names(tmod_res), names(cntr)))
    ret$tmod_res <- tmod_res
    ## XXX check that the tmod_res object is correct
  }

  message("checking annotation")
  if(!is.null(annot)) {
    ret$annot <- .check_df_primary_id(annot, primary_id)
  }

  message("returning")
  class(ret) <- c("seapiper_ds", class(ret))
  attr(ret, "primary_id") <- primary_id
  ret
}

#' S3 class for seapiper data sets
#'
#' S3 class for seapiper data sets
#' @return print.seapiper_ds does not return anything
#' @rdname seapiper_ds
#' @export
print.seapiper_ds <- function(x, ...) {
  .catf("Object of class seapiper_ds\n")
}



## make sure everything is where it is needed
.prepare_data <- function(x, primary_id, annot_default=NULL, tmod_dbs_default=NULL, tmod_map_default=NULL, save_memory=FALSE) {

  data <- imap(x, ~ {
    .id  <- .y
    #.prepare_data_single_pipeline(.id, .pip, primary_id, annot, cntr, tmod_res, tmod_dbs, save_memory=save_memory)
    if(is.null(.x$annot))    .x$annot    <- annot_default
    if(is.null(.x$tmod_map)) .x$tmod_map <- tmod_map_default
    if(is.null(.x$tmod_dbs)) .x$tmod_dbs <- tmod_dbs_default
    .x
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




