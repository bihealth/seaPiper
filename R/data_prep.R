## prepares the data structure for a single pipeline
.prepare_data_single_pipeline <- function(.id, .pip, primary_id, annot, cntr, tmod_res, tmod_dbs,
                                          contrast_cols=c(primary_id, "log2FoldChange", "pvalue", "padj")) {
  ret <- list()

  if(is.null(annot[[.id]])) {
    message(sprintf(" * Loading annotation for %s (consider using the annot option to speed this up)", .id))
    ret[["annot"]] <- get_annot(.pip)
  } else {
    ret[["annot"]] <- annot[[.id]]
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
  if(is.null(tmod_res[[.id]])) {
    message(sprintf(" * Loading tmod results for %s (consider using the tmod_res option to speed this up)", .id))
    ret[["tmod_res"]] <- tryCatch(
      get_tmod_res(.pip),
      error=function(e) {
        message(sprintf(" * tmod results unavailable for %s: %s", .id, conditionMessage(e)))
        NULL
      }
    )
  } else {
    ret[["tmod_res"]] <- tmod_res[[.id]]
  }

  if(is.null(tmod_dbs[[.id]])) {
    message(sprintf(" * Loading tmod databases for %s (consider using the tmod_dbs option to speed this up)", .id))
    ret[["tmod_dbs"]] <- tryCatch(
      get_tmod_dbs(.pip),
      error=function(e) {
        message(sprintf(" * tmod databases unavailable for %s: %s", .id, conditionMessage(e)))
        NULL
      }
    )
  } else {
    ret[["tmod_dbs"]] <- tmod_dbs[[.id]]
  }

 ## no longer necessary with the new tmod mset structure
 # ## get rid of unnecessary data
 #for(i in 1:length(ret[["tmod_dbs"]])) {
 #  ret[["tmod_dbs"]][[i]][["dbobj"]][["GENES2MODULES"]] <- NULL
 #}

  ## we only want the tmod object
  if(!is.null(ret[["tmod_dbs"]])) {
    ret[["tmod_dbs"]] <- map(ret[["tmod_dbs"]], ~ .x$dbobj)
  }

  ret[["tmod_map"]] <- tryCatch(
    get_tmod_mapping(.pip),
    error=function(e) {
      message(sprintf(" * tmod mapping unavailable for %s: %s", .id, conditionMessage(e)))
      NULL
    }
  )
  ret[["tmod_gl"]]  <- tryCatch(
    get_object(.pip, step="tmod", extension="gl.rds", as_list=TRUE),
    error=function(e) {
      message(sprintf(" * tmod gene lists unavailable for %s: %s", .id, conditionMessage(e)))
      NULL
    }
  )

  ## we only need the order of the genes from one tmod db
  if(!is.null(ret[["tmod_gl"]]) && !is.null(ret[["annot"]]) && primary_id %in% colnames(ret[["annot"]])) {
    ret[["tmod_gl"]] <- map(ret[["tmod_gl"]], ~  # one for each of contrast
                               map(.[[1]], ~  # one for each sorting type
                                   match(names(.), ret[["annot"]][[primary_id]])))
  }

  ret[["config"]]   <- get_config(.pip)
  ret[["covar"]]    <- get_covariates(.pip)

  if(is.null(ret[["config"]][["dataset_title"]])) {
    ret[["config"]][["dataset_title"]] <- .id
  }

  if(!is.null(ret[["annot"]])) {
    ret[["annot_linkout"]] <- .prep_annot_linkout(ret[["annot"]], ret[["config"]])
  } else {
    ret[["annot_linkout"]] <- NULL
  }

  if(!is.null(ret[["tmod_dbs"]])) {
    ret[["dbs"]]     <- names(ret[["tmod_dbs"]])
    ret[["sorting"]] <- ret[["config"]]$tmod$sort_by
  } else {
    ret[["dbs"]]     <- NULL
    ret[["sorting"]] <- NULL
  }

  ret[["rld"]] <- get_object(.pip, step="DESeq2", extension="rld.blind.rds")
  ret[["rld"]] <- .extract_rld_matrix(ret[["rld"]], .id)

  ## prepare the PCA
  mtx <- t(ret[["rld"]])
  vars <- apply(mtx, 2, var)
  sel  <- vars > 1e-26
  ret[["pca"]] <- prcomp(mtx[ , sel], scale.=TRUE)$x
  
  if(!is.null(ret[["cntr"]]) && !is.null(ret[["config"]]$contrasts$contrast_list)) {
    ret[["cntr_titles"]]        <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "ID")
    names(ret[["cntr_titles"]]) <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "title")
    ret[["cntr_titles"]]        <- ret[["cntr_titles"]][ ret[["cntr_titles"]] %in% names(ret[["cntr"]]) ]
  } else {
    ret[["cntr_titles"]] <- NULL
  }

  ret
}

# Extract an expression matrix from rld-like objects without Bioconductor generics.
# Returns a numeric matrix (or stops with an informative error).
.extract_rld_matrix <- function(x, dataset_id) {
  if(is.null(x)) {
    stop(sprintf("rld object is missing for dataset `%s`", dataset_id))
  }

  if(is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if(is.matrix(x)) {
    return(x)
  }

  if(methods::isS4(x) && "assays" %in% methods::slotNames(x)) {
    assays_obj <- methods::slot(x, "assays")
    if(methods::isS4(assays_obj) && "data" %in% methods::slotNames(assays_obj)) {
      assays_obj <- methods::slot(assays_obj, "data")
    }
    assays_list <- tryCatch(as.list(assays_obj), error=function(e) NULL)
    if(!is.null(assays_list) && length(assays_list) > 0) {
      first <- assays_list[[1]]
      if(is.data.frame(first)) {
        first <- as.matrix(first)
      }
      if(is.matrix(first)) {
        return(first)
      }
    }
  }

  stop(
    sprintf(
      "Cannot extract matrix from rld object for dataset `%s`; expected matrix/data.frame or S4 object with assay data",
      dataset_id
    )
  )
}


## build linkout templates for available annotation columns
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
