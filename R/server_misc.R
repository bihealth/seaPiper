## resolve a single sample-id column from stored dataset metadata
.resolve_sample_id_col <- function(sample_id, default="SampleID") {
  if(is.null(sample_id)) {
    return(default)
  }

  if(is.list(sample_id)) {
    sample_id <- unlist(sample_id, use.names=FALSE)
  }

  sample_id <- as.character(sample_id)
  sample_id <- sample_id[!is.na(sample_id) & nzchar(sample_id)]
  sample_id <- unique(sample_id)

  if(length(sample_id) == 0L) {
    return(default)
  }

  if(length(sample_id) > 1L) {
    stop("All datasets must use the same `sample_id` value.")
  }

  sample_id[[1L]]
}

## normalize covariates for expression modules and enforce sample-id consistency
.normalize_expression_covar <- function(covar, sample_id_col="SampleID") {
  if(is.null(covar)) {
    return(NULL)
  }

  covar_list <- covar
  if(is.data.frame(covar_list)) {
    covar_list <- list(default=covar_list)
  }

  if(!is.list(covar_list) || length(covar_list) == 0L) {
    stop("`covar` must be a data frame or a non-empty list of data frames.")
  }

  if(is.null(names(covar_list))) {
    names(covar_list) <- paste0("dataset_", seq_along(covar_list))
  }

  out <- covar_list
  for(i in seq_along(covar_list)) {
    ds <- names(covar_list)[i]
    covar_df <- covar_list[[i]]

    if(!is.data.frame(covar_df)) {
      stop(sprintf("Covariates for dataset `%s` must be a data frame.", ds))
    }

    if(!sample_id_col %in% colnames(covar_df)) {
      stop(sprintf(
        "Dataset `%s`: covariates must contain sample ID column `%s`.",
        ds,
        sample_id_col
      ))
    }

    sample_ids <- as.character(covar_df[[sample_id_col]])
    if(any(is.na(sample_ids) | sample_ids == "")) {
      stop(sprintf(
        "Dataset `%s`: `%s` contains missing or empty sample IDs.",
        ds,
        sample_id_col
      ))
    }

    if(anyDuplicated(sample_ids)) {
      stop(sprintf(
        "Dataset `%s`: `%s` must contain unique sample IDs.",
        ds,
        sample_id_col
      ))
    }

    rownames(covar_df) <- sample_ids
    out[[i]] <- covar_df
  }

  out
}

## server logic for heatmap, disco, volcano, and pca modules
.seapiper_server_misc <- function(input, output, session, data, gene_id,
                                  selection=NULL, enable_disco=FALSE,
                                  enable_volcano=FALSE, enable_heatmap=FALSE,
                                  enable_pca=FALSE, palettes=NULL, covar=NULL,
                                  sample_id_col="SampleID") {
  if(is.null(covar)) {
    covar <- data[["covar"]]
  }

  if(isTRUE(enable_heatmap)) {
    bioshmods::heatmapServer("heatmap",
                             annot=data[["annot"]],
                             exprs=data[["rld"]],
                             cntr=data[["cntr"]],
                             covar=covar,
                             selection=selection,
                             palettes=palettes,
                             sample_id_col=sample_id_col)
  }

  if(isTRUE(enable_disco)) {
    discoServer("disco", data[["cntr"]], data[["annot"]], gene_id=gene_id)
  }

  if(isTRUE(enable_volcano)) {
    volcanoServer("volcano", data[["cntr"]], annot=data[["annot"]],
                  ui_config = list(show_button_label = "Show heatmap"),
                  gene_id=gene_id, selection=selection)
  }

  if(isTRUE(enable_pca)) {
    pcaServer("pca", data[["pca"]], covar)
  }
}
