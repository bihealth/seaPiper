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
seapiperdata_from_objects <- function(cntr, annot, exprs,
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
    tmod_res = NULL,
    tmod_dbs = NULL,
    tmod_map = NULL,
    tmod_gl  = NULL,
    config   = list(default=config),
    covar    = list(default=covar),
    annot_linkout = list(default=.prep_annot_linkout(annot, config)),
    dbs      = NULL,
    sorting  = NULL,
    rld      = list(default=exprs),
    pca      = list(default=pca),
    cntr_titles = list(default=cntr_titles)
  )

  .as_seapiperdata(data)
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
seapiperdata_from_rseasnap <- function(pip, primary_id="PrimaryID",
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

  .as_seapiperdata(transpose(data))
}


#' Build seaPiperData from a YAML file of RDS paths
#'
#' Load a `seaPiperData` object from a YAML file that points to serialized R
#' objects (`.rds`). Relative paths are resolved relative to the YAML file
#' itself.
#'
#' The YAML must define one or more datasets in a top-level `datasets`
#' mapping.
#'
#' Supported dataset keys and expected `.rds` object types are:
#' * `annot`: `data.frame` with gene annotation (optionally includes `primary_id`)
#' * `cntr`: named `list` of contrast `data.frame`s
#' * `tmod_res`: named `list` of tmod results (per contrast)
#' * `tmod_dbs`: named `list` of tmod DB objects (or wrappers with `dbobj`)
#' * `tmod_map`: tmod mapping object (`get_tmod_mapping` output)
#' * `tmod_gl`: tmod gene-ranking object (`gl.rds`-like nested structure)
#' * `config`: `list` with sea-snap style config fields
#' * `covar`: `data.frame` of sample covariates (required)
#' * `rld`: numeric matrix/data.frame (genes x samples)
#' * `pca`: PCA scores matrix/data.frame (samples x PCs)
#'
#' @param yaml_file path to a YAML file describing RDS locations
#' @param primary_id name of the primary gene ID column
#' @examples
#' # # Example YAML structure (all paths are relative to this YAML file):
#' # datasets:
#' #   default:
#' #     annot: annot.rds
#' #     cntr: cntr.rds
#' #     tmod_res: tmod_res.rds
#' #     tmod_dbs: tmod_dbs.rds
#' #     tmod_map: tmod_map.rds
#' #     tmod_gl: tmod_gl.rds
#' #     config: config.rds
#' #     covar: covar.rds
#' #     rld: rld.rds
#' @importFrom yaml read_yaml
#' @export
seapiperdata_from_yaml <- function(yaml_file, primary_id="PrimaryID") {
  stopifnot(is.character(yaml_file), length(yaml_file) == 1)
  if(!file.exists(yaml_file)) {
    stop(sprintf("YAML file not found: %s", yaml_file))
  }

  spec <- read_yaml(yaml_file)
  if(is.null(spec) || !is.list(spec)) {
    stop("`yaml_file` must parse to a YAML mapping")
  }

  datasets <- spec$datasets %||% list()
  if(is.null(datasets) || !is.list(datasets) || length(datasets) == 0) {
    stop("YAML must define a top-level `datasets` mapping with at least one dataset")
  }

  if(is.null(names(datasets)) || any(names(datasets) == "")) {
    stop("All entries under `datasets` must be named by dataset ID")
  }

  base_dir <- dirname(normalizePath(yaml_file, mustWork=TRUE))
  data <- imap(datasets, ~ .load_custom_yaml_dataset(.x, .y, base_dir, primary_id))

  .as_seapiperdata(transpose(data))
}

#' Merge two seaPiperData objects
#'
#' Combine two `seaPiperData` objects into one. Dataset IDs from both inputs are
#' preserved. If the same dataset ID appears in both objects, IDs are renamed
#' using `suffixes` (similarly to `merge()`), e.g. `default.x` and `default.y`.
#'
#' @param x first `seaPiperData` object
#' @param y second `seaPiperData` object
#' @param suffixes character vector of length 2 used when a dataset ID appears
#'   in both objects.
#' @return A merged object of class `seaPiperData`.
#' @export
merge_seapiperdata <- function(x, y, suffixes=c(".x", ".y")) {
  .assert_seapiperdata(x, "x")
  .assert_seapiperdata(y, "y")

  if(!is.character(suffixes) || length(suffixes) != 2 || anyNA(suffixes)) {
    stop("`suffixes` must be a character vector of length 2")
  }

  ids_x <- .dataset_ids_from_seapiperdata(x)
  ids_y <- .dataset_ids_from_seapiperdata(y)
  renaming <- .dataset_renaming(ids_x, ids_y, suffixes)
  fields <- .merge_field_names(x, y)

  merged <- map(fields, ~ {
    x_values <- .normalize_field_values(x[[.x]], ids_x, .x)
    y_values <- .normalize_field_values(y[[.x]], ids_y, .x)
    x_values <- .rename_field_values(x_values, renaming$map_x)
    y_values <- .rename_field_values(y_values, renaming$map_y)
    combined <- c(x_values, y_values)
    if(length(combined) > 0 && all(vapply(combined, is.null, logical(1)))) {
      NULL
    } else {
      combined
    }
  })
  names(merged) <- fields

  .as_seapiperdata(merged)
}

#' Build seaPiperData from basic inputs
#'
#' Deprecated alias of `seapiperdata_from_objects()`.
#' @param ... arguments passed to `seapiperdata_from_objects()`
#' @export
make_seapiperdata <- function(...) {
  .Deprecated("seapiperdata_from_objects")
  seapiperdata_from_objects(...)
}

#' Convert Rseasnap pipeline output into seaPiperData
#'
#' Deprecated alias of `seapiperdata_from_rseasnap()`.
#' @param ... arguments passed to `seapiperdata_from_rseasnap()`
#' @export
rseasnap_to_seapiperdata <- function(...) {
  .Deprecated("seapiperdata_from_rseasnap")
  seapiperdata_from_rseasnap(...)
}

#' Build seaPiperData from a YAML file of RDS paths
#'
#' Deprecated alias of `seapiperdata_from_yaml()`.
#' @param ... arguments passed to `seapiperdata_from_yaml()`
#' @export
custom_yaml_to_seapiperdata <- function(...) {
  .Deprecated("seapiperdata_from_yaml")
  seapiperdata_from_yaml(...)
}

# Validate that an input is a seaPiperData object.
# Returns NULL; used for side effects (errors).
.assert_seapiperdata <- function(x, arg) {
  if(!inherits(x, "seaPiperData")) {
    stop(sprintf("`%s` must be a seaPiperData object", arg))
  }
}

# Compute the union of dataset IDs from one or more seaPiperData objects.
# Returns a character vector of dataset IDs (possibly empty).
.dataset_ids_from_seapiperdata <- function(...) {
  objects <- list(...)
  ids <- character(0)

  for(obj in objects) {
    for(field in names(obj)) {
      value <- obj[[field]]
      if(is.list(value) && !is.null(names(value)) && length(value) > 0) {
        ids <- c(ids, names(value))
      }
    }
  }

  unique(ids)
}

# Determine merged top-level field names in stable schema-first order.
# Returns character vector of field names to merge.
.merge_field_names <- function(x, y) {
  schema_fields <- c("annot", "cntr", "tmod_res", "tmod_dbs", "tmod_map", "tmod_gl",
                     "config", "covar", "annot_linkout", "dbs", "sorting",
                     "rld", "pca", "cntr_titles")
  unique(c(schema_fields, names(x), names(y)))
}

# Build dataset-ID renaming maps for two objects using suffixes for overlaps.
# Returns a list with `map_x` and `map_y` named character vectors.
.dataset_renaming <- function(ids_x, ids_y, suffixes) {
  overlap <- intersect(ids_x, ids_y)

  map_x <- setNames(ids_x, ids_x)
  map_y <- setNames(ids_y, ids_y)

  map_x[overlap] <- paste0(overlap, suffixes[[1]])
  map_y[overlap] <- paste0(overlap, suffixes[[2]])

  .assert_unique_target_ids(unname(map_x), "x")
  .assert_unique_target_ids(unname(map_y), "y")
  .assert_unique_target_ids(c(unname(map_x), unname(map_y)), "merged result")

  list(map_x=map_x, map_y=map_y)
}

# Validate that renamed dataset IDs are unique.
# Returns NULL; used for side effects (errors).
.assert_unique_target_ids <- function(ids, label) {
  dup <- unique(ids[duplicated(ids)])
  if(length(dup) > 0) {
    stop(sprintf(
      "Duplicate dataset IDs in %s after renaming: %s",
      label,
      paste(dup, collapse=", ")
    ))
  }
}

# Normalize one top-level field to a named list indexed by dataset IDs.
# Returns a list of length `length(dataset_ids)` with names `dataset_ids`.
.normalize_field_values <- function(values, dataset_ids, field_name) {
  ret <- setNames(vector("list", length(dataset_ids)), dataset_ids)

  if(is.null(values)) {
    return(ret)
  }
  if(!is.list(values)) {
    stop(sprintf("Field `%s` must be a list or NULL", field_name))
  }
  if(is.null(names(values)) || any(names(values) == "")) {
    stop(sprintf("Field `%s` must be a named list", field_name))
  }

  keep <- names(values) %in% dataset_ids
  ret[names(values)[keep]] <- values[keep]
  ret
}

# Rename one normalized field list according to a dataset-ID map.
# Returns a named list with updated dataset IDs.
.rename_field_values <- function(values, id_map) {
  if(length(values) == 0) {
    return(values)
  }
  names(values) <- unname(id_map[names(values)])
  values
}

# Load and normalize one dataset section from the YAML spec.
# Returns a single normalized dataset list matching the seaPiperData schema.
.load_custom_yaml_dataset <- function(dataset_spec, dataset_name, base_dir, primary_id) {
  if(!is.list(dataset_spec)) {
    stop(sprintf("Dataset `%s` must be a YAML mapping", dataset_name))
  }

  keys <- c("annot", "cntr", "tmod_res", "tmod_dbs", "tmod_map", "tmod_gl",
            "config", "covar", "rld", "pca")
  ret <- .read_dataset_keys(dataset_spec, keys, base_dir, dataset_name)
  .validate_required_and_report_missing(ret, keys, dataset_name)
  ret <- .normalize_core_inputs(ret, dataset_name, primary_id)
  ret <- .ensure_dataset_config(ret, dataset_name)
  ret <- .normalize_covar(ret)
  ret <- .normalize_tmod_inputs(ret, dataset_name, primary_id)
  ret <- .compute_pca_if_missing(ret, dataset_name)
  .finalize_dataset(ret)
}

# Load configured keys for one dataset.
# Returns a list with one element per supported dataset key.
.read_dataset_keys <- function(dataset_spec, keys, base_dir, dataset_name) {
  ret <- list()
  for(key in keys) {
    ret[[key]] <- .yaml_spec_value(dataset_spec[[key]], base_dir, dataset_name, key)
  }
  ret
}

# Validate required keys and report which optional keys were omitted.
# Returns NULL; used for side effects (errors/messages).
.validate_required_and_report_missing <- function(ret, keys, dataset_name) {
  if(is.null(ret$covar)) {
    stop(sprintf("Dataset `%s`: `covar` is required", dataset_name))
  }

  optional_keys <- keys[keys != "covar"]
  missing_optional <- optional_keys[vapply(ret[optional_keys], is.null, logical(1))]
  if(length(missing_optional) > 0) {
    message(
      sprintf(
        "Dataset `%s`: optional keys not provided: %s. Related app modules may be disabled.",
        dataset_name,
        paste(missing_optional, collapse=", ")
      )
    )
  }
}

# Normalize annot/cntr/rld inputs and validate structural requirements.
# Returns updated dataset list with normalized core inputs.
.normalize_core_inputs <- function(ret, dataset_name, primary_id) {
  if(!is.null(ret$annot)) {
    ret$annot <- as.data.frame(ret$annot)
    if(!primary_id %in% colnames(ret$annot)) {
      if(!is.null(rownames(ret$annot))) {
        ret$annot[[primary_id]] <- rownames(ret$annot)
      } else {
        stop(sprintf("Dataset `%s`: `annot` needs `%s` column or row names", dataset_name, primary_id))
      }
    }
  }

  if(!is.null(ret$cntr)) {
    if(!is.list(ret$cntr)) {
      stop(sprintf("Dataset `%s`: `cntr` must be a named list", dataset_name))
    }
    if(is.null(names(ret$cntr))) {
      names(ret$cntr) <- paste0("contrast_", seq_along(ret$cntr))
    }
    ret$cntr <- map(ret$cntr, ~ {
      .x <- as.data.frame(.x)
      if(!primary_id %in% colnames(.x)) {
        if(!is.null(rownames(.x))) {
          .x[[primary_id]] <- rownames(.x)
        } else {
          stop(sprintf("Dataset `%s`: each contrast needs `%s` column or row names", dataset_name, primary_id))
        }
      }
      .x
    })
  }

  if(!is.null(ret$rld)) {
    if(is.data.frame(ret$rld)) {
      ret$rld <- as.matrix(ret$rld)
    }
    if(!is.matrix(ret$rld)) {
      stop(sprintf("Dataset `%s`: `rld` must be a matrix/data.frame", dataset_name))
    }
    if(is.null(colnames(ret$rld))) {
      stop(sprintf("Dataset `%s`: `rld` must have column names", dataset_name))
    }
    if(is.null(rownames(ret$rld))) {
      stop(sprintf("Dataset `%s`: `rld` must have row names", dataset_name))
    }
  }

  ret
}

# Ensure config is present and carries dataset_title.
# Returns dataset list with a usable config object.
.ensure_dataset_config <- function(ret, dataset_name) {
  if(is.null(ret$config)) {
    ret$config <- list(
      organism   = list(name="unknown", taxon="unknown"),
      experiment = list(design_formula="not_set"),
      contrasts  = list(contrast_list=map(names(ret$cntr %||% list()), ~ list(ID=.x, title=.x))),
      tmod       = list(databases=list(), sort_by=character(0)),
      filter     = list(low_counts=NA_integer_, min_counts=NA_integer_, min_count_n=NA_integer_),
      dataset_title = dataset_name
    )
    message(sprintf("Dataset `%s`: `config` not provided; using default config.", dataset_name))
  } else if(is.null(ret$config$dataset_title)) {
    ret$config$dataset_title <- dataset_name
  }

  ret
}

# Normalize covariate table and row names.
# Returns dataset list with normalized covar data.frame.
.normalize_covar <- function(ret) {
  ret$covar <- as.data.frame(ret$covar)
  if(is.null(rownames(ret$covar)) && "ID" %in% colnames(ret$covar)) {
    rownames(ret$covar) <- ret$covar[["ID"]]
  }
  ret
}

# Normalize tmod-related structures and gene-list index mapping.
# Returns dataset list with normalized tmod_dbs/tmod_gl.
.normalize_tmod_inputs <- function(ret, dataset_name, primary_id) {
  if(!is.null(ret$tmod_dbs) && is.list(ret$tmod_dbs) && length(ret$tmod_dbs) > 0) {
    has_dbobj <- vapply(ret$tmod_dbs, function(.x) {
      is.list(.x) && !is.null(.x$dbobj)
    }, logical(1))
    if(all(has_dbobj)) {
      ret$tmod_dbs <- map(ret$tmod_dbs, ~ .x$dbobj)
    }
  }

  if(!is.null(ret$tmod_gl) && !is.null(ret$annot) && primary_id %in% colnames(ret$annot)) {
    ret$tmod_gl <- .convert_tmod_gl(ret$tmod_gl, ret$annot[[primary_id]])
  } else if(!is.null(ret$tmod_gl) && is.null(ret$annot)) {
    message(sprintf("Dataset `%s`: `tmod_gl` provided without `annot`; loaded as-is.", dataset_name))
  }

  ret
}

# Compute PCA from rld when not explicitly provided.
# Returns dataset list with pca populated if computable.
.compute_pca_if_missing <- function(ret, dataset_name) {
  if(is.null(ret$pca) && !is.null(ret$rld)) {
    mtx  <- t(ret$rld)
    vars <- apply(mtx, 2, var)
    sel  <- vars > 1e-26
    ret$pca <- prcomp(mtx[ , sel], scale.=TRUE)$x
    message(sprintf("Dataset `%s`: computed `pca` from `rld`.", dataset_name))
  }
  ret
}

# Fill derived fields and finalize default fallbacks for seaPiperData schema.
# Returns a fully normalized dataset entry.
.finalize_dataset <- function(ret) {
  if(!is.null(ret$annot)) {
    ret$annot_linkout <- .prep_annot_linkout(ret$annot, ret$config)
  } else {
    ret$annot_linkout <- NULL
  }
  if(!is.null(ret$tmod_dbs)) {
    ret$dbs <- names(ret$tmod_dbs)
    ret$sorting <- ret$config$tmod$sort_by %||% NULL
  } else {
    ret$dbs <- NULL
    ret$sorting <- NULL
  }

  if(!is.null(ret$cntr) && !is.null(ret$config$contrasts$contrast_list)) {
    ret$cntr_titles <- map_chr(ret$config$contrasts$contrast_list, `[[`, "ID")
    names(ret$cntr_titles) <- map_chr(ret$config$contrasts$contrast_list, `[[`, "title")
    ret$cntr_titles <- ret$cntr_titles[ ret$cntr_titles %in% names(ret$cntr) ]
  } else if(!is.null(ret$cntr)) {
    ret$cntr_titles <- setNames(names(ret$cntr), names(ret$cntr))
  } else {
    ret$cntr_titles <- NULL
  }

  ret
}

# Normalize constructor output to seaPiperData and collapse globally-absent sections.
# Returns a seaPiperData object where fields that are NULL across all datasets are top-level NULL.
.as_seapiperdata <- function(data) {
  data <- .collapse_all_null_fields(data)
  structure(data, class="seaPiperData")
}

# Collapse list fields for which every dataset entry is NULL.
# Returns the input list with all-NULL fields replaced by top-level NULL.
.collapse_all_null_fields <- function(data) {
  for(name in names(data)) {
    value <- data[[name]]
    if(is.list(value) && length(value) > 0 && all(vapply(value, is.null, logical(1)))) {
      data[[name]] <- NULL
    }
  }
  data
}

# Resolve one YAML key value: load from RDS path or pass through literal object.
# Returns either NULL, a deserialized R object, or the original non-path value.
.yaml_spec_value <- function(value, base_dir, dataset_name, key) {
  if(is.null(value)) {
    return(NULL)
  }

  if(is.character(value) && length(value) == 1) {
    path <- .resolve_path(value, base_dir)
    if(!file.exists(path)) {
      stop(sprintf("Dataset `%s`: file for `%s` not found: %s", dataset_name, key, path))
    }
    return(readRDS(path))
  }

  value
}

# Null-coalescing helper for concise default handling.
# Returns left-hand side unless it is NULL, otherwise returns right-hand side.
`%||%` <- function(lhs, rhs) {
  if(is.null(lhs)) rhs else lhs
}

# Resolve a potentially relative file path against the YAML directory.
# Returns a normalized path string.
.resolve_path <- function(path, base_dir) {
  if(.is_absolute_path(path)) {
    return(path.expand(path))
  }
  normalizePath(file.path(base_dir, path), mustWork=FALSE)
}

# Test whether a path is absolute (Unix, home, or Windows drive path).
# Returns TRUE/FALSE.
.is_absolute_path <- function(path) {
  grepl("^(/|~|[A-Za-z]:[\\\\/])", path)
}

# Convert tmod gene-list objects from named genes to annotation index vectors.
# Returns a recursively transformed object with matched integer indices.
.convert_tmod_gl <- function(x, annot_ids) {
  if(is.list(x)) {
    return(map(x, ~ .convert_tmod_gl(.x, annot_ids)))
  }

  nms <- names(x)
  if(is.null(nms) || length(nms) == 0) {
    return(x)
  }

  match(nms, annot_ids)
}
