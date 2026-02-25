# Coerce one object to data.frame for tabular export.
# Returns a data.frame or NULL when the object is not coercible.
.as_export_df <- function(x) {
  if(is.data.frame(x)) {
    return(x)
  }
  if(is.matrix(x)) {
    return(as.data.frame(x, stringsAsFactors=FALSE))
  }
  NULL
}

# Normalize one seaPiper field to fileExport-compatible dataframe payloads.
# Returns a data.frame, list of data.frames, list of lists of data.frames, or NULL.
.as_export_dataframe_payload <- function(x) {
  if(is.null(x)) {
    return(NULL)
  }

  top_df <- .as_export_df(x)
  if(!is.null(top_df)) {
    return(top_df)
  }

  if(!is.list(x) || length(x) == 0L) {
    return(NULL)
  }

  lvl1 <- lapply(x, .as_export_df)
  keep_lvl1 <- !vapply(lvl1, is.null, logical(1))
  if(any(keep_lvl1)) {
    lvl1 <- lvl1[keep_lvl1]
    if(length(lvl1) > 0L) {
      return(lvl1)
    }
  }

  lvl2 <- lapply(x, function(.x) {
    if(!is.list(.x) || length(.x) == 0L) {
      return(NULL)
    }

    .inner <- lapply(.x, .as_export_df)
    keep_inner <- !vapply(.inner, is.null, logical(1))
    .inner <- .inner[keep_inner]
    if(length(.inner) == 0L) {
      return(NULL)
    }
    .inner
  })

  keep_lvl2 <- !vapply(lvl2, is.null, logical(1))
  lvl2 <- lvl2[keep_lvl2]
  if(length(lvl2) > 0L) {
    return(lvl2)
  }

  NULL
}

# Append one export specification if data are present.
# Returns updated list of export specification entries.
.append_export_spec <- function(specs, data, title, description, data_type, help=NULL) {
  is_empty_list <- is.list(data) && length(data) == 0L
  if(is.null(data) || is_empty_list) {
    return(specs)
  }

  specs[[length(specs) + 1L]] <- list(
    data=data,
    title=title,
    description=description,
    data_type=data_type,
    help=help
  )

  specs
}


## recursively collapse a list of dataframes
.collapse_nested_ldf <- function(nested_ldf, col_names) {

  if(is.null(nested_ldf)) {
    return(NULL)
  }

  if(is.data.frame(nested_ldf)) {
    return(nested_ldf)
  }

  if(length(col_names) == 0L) {
    stop("No more column names provided for collapsing nested list.")
  }

  cname <- col_names[1L]
  col_names <- col_names[-1L]

  if(!is.list(nested_ldf) || length(nested_ldf) == 0L) {
    stop("Expected a non-empty list for collapsing, but got an empty list or non-list object.")
  }

  res <- lapply(names(nested_ldf), function(name) {
    inner <- nested_ldf[[name]]
    if(is.list(inner)) {
      inner <- .collapse_nested_ldf(inner, col_names=col_names)
    }
    if(!is.data.frame(inner)) {
      stop("Expected a data frame at the innermost level, but got a non-data frame object.")
    }

    inner[[cname]] <- name
    inner <- inner[, c(ncol(inner), 1:(ncol(inner) - 1)), drop=FALSE]
  })

  res <- do.call(rbind, res)
  res
}


# Coerce tmod result lists to a list of lists of data frames for export.
.tmod_res_as_lldf <- function(tmod_res) {
  if(is.null(tmod_res)) {
    return(NULL)
  }

  if(!is.list(tmod_res) || length(tmod_res) == 0L) {
    return(NULL)
  }

  print(names(tmod_res))
  res <- lapply(tmod_res, function(tmod_single_dataset) {
                  lapply(tmod_single_dataset, function(tmod_single_cntr) {
                    .collapse_nested_ldf(tmod_single_cntr, c("Database", "Sorting"))
                  })
  })

  # this should now be a list of data frames
  res

}

# Build bioshmods::fileExport object specifications from seaPiperData.
# Returns a non-empty list ready for fileExportUI/fileExportServer.
.build_export_objects <- function(data) {
  specs <- list()

  specs <- .append_export_spec(
    specs,
    .as_export_dataframe_payload(data[["cntr"]]),
    title="Differential Expression Contrasts",
    description="Per-dataset contrast tables (xlsx/zip export)",
    data_type="dataframes",
    help="Main gene-level differential expression results shown across Gene browser, Volcano, and Disco."
  )

  specs <- .append_export_spec(
    specs,
    .as_export_dataframe_payload(data[["annot"]]),
    title="Gene Annotation",
    description="Per-dataset annotation tables",
    data_type="dataframes",
    help="Gene annotation used for labels, lookups, and link-outs."
  )

  specs <- .append_export_spec(
    specs,
    .as_export_dataframe_payload(data[["covar"]]),
    title="Covariates",
    description="Per-dataset sample metadata",
    data_type="dataframes",
    help="Sample-level covariates shown in Workflow Info and used by plotting modules."
  )

  specs <- .append_export_spec(
    specs,
    .as_export_dataframe_payload(data[["rld"]]),
    title="Expression Matrix (rlog)",
    description="Per-dataset expression matrices (genes x samples)",
    data_type="dataframes",
    help="Regularized log expression values used in Gene browser plotting."
  )

  specs <- .append_export_spec(
    specs,
    .as_export_dataframe_payload(data[["pca"]]),
    title="PCA Coordinates",
    description="Per-dataset PCA scores",
    data_type="dataframes",
    help="PCA coordinates displayed in the PCA module."
  )

  specs <- .append_export_spec(
    specs,
    .tmod_res_as_lldf(data[["tmod_res"]]),
    title="Tmod Results",
    description="Gene set enrichment results as Excel tables",
    data_type="dataframes",
    help="Tmod enrichment results for all contrasts and databases, pre-processed into a list of lists of data frames for export. Each top-level list corresponds to a contrast, and contains a list of data frames for each database, with an additional column indicating the sorting method (e.g. 'pval')."
  )

  specs <- .append_export_spec(
    specs,
    data[["tmod_dbs"]],
    title="Tmod Databases (Raw Object)",
    description="Gene set database objects as RDS",
    data_type="rds",
    help="Raw tmod database objects used for enrichment statistics and module details."
  )

  specs <- .append_export_spec(
    specs,
    data[["tmod_map"]],
    title="Tmod Mapping (Raw Object)",
    description="Module-to-gene mapping as RDS",
    data_type="rds",
    help="Raw tmod mapping object."
  )

  specs <- .append_export_spec(
    specs,
    data[["tmod_gl"]],
    title="Tmod Gene Lists (Raw Object)",
    description="Precomputed ranked gene lists as RDS",
    data_type="rds",
    help="Raw tmod gene ranking object used for evidence plots."
  )

  specs <- .append_export_spec(
    specs,
    data[["config"]],
    title="Pipeline Configuration (Raw Object)",
    description="Per-dataset configuration as RDS",
    data_type="rds",
    help="Contains experiment design, filters, and contrast metadata."
  )

  specs <- .append_export_spec(
    specs,
    data,
    title="Complete seaPiperData Object",
    description="Full application data object as RDS",
    data_type="rds",
    help="Full object backup preserving all currently loaded structures."
  )

  specs
}
