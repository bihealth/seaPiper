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
#' @param tmod_res Results of the gene set enrichment analysis. If the list is unnamed, then
#' names will be generated (ID1, ID2, ...). The sets of names of `tmod_res`
#' and `cntr` must be equal. Each element of the list must
#' itself be either a data frame with gene set enrichment results or a named list of data frames
#' with gene set enrichment results
#' @param annot Annotation data frame. If defined, the column indicated by
#' the `primary_id` parameter must be present; alternatively, the rows
#' define the primary gene IDs.
#' @param primary_id Name of the column in contrast data frames which holds
#' the primary gene identifier.
#' @return An object of the class seapiper_ds
#' @export
# XXX how should we deal with tmod gene sets results? biglist or list of
# dataframes?
new_seapiper_dataset <- function(
                                 cntr,
                                 tmod_res=NULL,
                                 annot=NULL,
                                 primary_id = "PrimaryID"

  
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

  cntr <- imap(cntr, ~ {
    # XXX check for .x being a DF
    .x <- .check_df_primary_id(.x, primary_id)
    .x
  })

  ret$cntr <- cntr

  if(!is.null(tmod_res)) {
    if(is.null(names(tmod_res))) {
      names(tmod_res) <- paste0("ID", seq_along(tmod_res))
    }

    stopifnot(setequal(names(tmod_res), names(cntr)))
    ret$tmod_res <- tmod_res
    ## XXX check that the tmod_res object is correct
  }

  if(!is.null(annot)) {
    ret$annot <- .check_df_primary_id(annot, primary_id)
  }

  class(ret) <- c("seapiper_ds", class(ret))
  attr(ret, "primary_id") <- primary_id
  ret
}

#' S3 class for seapiper data sets
#'
#' S3 class for seapiper data sets
#' @return print.seapiper_ds does not return anything
#' @rdname seapiper_ds-class
#' @export
print.seapiper_ds <- function(x, ...) {
  .catf("Object of class seapiper_ds\n")
}
