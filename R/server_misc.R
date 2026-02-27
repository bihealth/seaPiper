## detect a shared sample-id column suitable for bioshmods::heatmapServer
.detect_heatmap_sample_id_col <- function(covar) {
  candidates <- c("SampleID", "ID", "label", "sample", "sample_id")

  if(is.null(covar)) {
    return(NULL)
  }

  if(is.data.frame(covar)) {
    hit <- candidates[candidates %in% colnames(covar)]
    if(length(hit) > 0L) {
      return(hit[[1L]])
    }
    return(NULL)
  }

  if(is.list(covar) && length(covar) > 0L) {
    covar_df <- covar[vapply(covar, is.data.frame, logical(1))]
    if(length(covar_df) == 0L) {
      return(NULL)
    }

    for(cn in candidates) {
      if(all(vapply(covar_df, function(df) cn %in% colnames(df), logical(1)))) {
        return(cn)
      }
    }
  }

  NULL
}

## server logic for heatmap, disco, volcano, and pca modules
.seapiper_server_misc <- function(input, output, session, data, gene_id, enable_disco=FALSE,
                                  enable_volcano=FALSE, enable_heatmap=FALSE,
                                  enable_pca=FALSE, palettes=NULL) {
  if(isTRUE(enable_heatmap)) {
    sample_id_col <- .detect_heatmap_sample_id_col(data[["covar"]])
    covar <- data[["covar"]]
    if(is.null(sample_id_col)) {
      covar <- NULL
    }

    bioshmods::heatmapServer("heatmap",
                             annot=data[["annot"]],
                             exprs=data[["rld"]],
                             cntr=data[["cntr"]],
                             covar=covar,
                             palettes=palettes,
                             sample_id_col=if(is.null(sample_id_col)) "SampleID" else sample_id_col)
  }

  if(isTRUE(enable_disco)) {
    discoServer("disco", data[["cntr"]], data[["annot"]], gene_id=gene_id)
  }

  if(isTRUE(enable_volcano)) {
    volcanoServer("volcano", data[["cntr"]], annot=data[["annot"]], gene_id=gene_id)
  }

  if(isTRUE(enable_pca)) {
    pcaServer("pca", data[["pca"]], data[["covar"]])
  }
}
