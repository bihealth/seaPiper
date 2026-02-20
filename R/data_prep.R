## load missing objects, prepare the data etc.
.prepare_data <- function(pip, primary_id, annot=NULL, cntr=NULL, tmod_res=NULL, tmod_dbs=NULL) {

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

  return(transpose(data))
}

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

 ## no longer necessary with the new tmod mset structure
 # ## get rid of unnecessary data
 #for(i in 1:length(ret[["tmod_dbs"]])) {
 #  ret[["tmod_dbs"]][[i]][["dbobj"]][["GENES2MODULES"]] <- NULL
 #}

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

  ret[["dbs"]]     <- names(ret[["tmod_dbs"]])
  ret[["sorting"]] <- ret[["config"]]$tmod$sort_by

  ret[["rld"]]     <- get_object(.pip, step="DESeq2", extension="rld.blind.rds")
  ret[["rld"]]     <- assay(ret[["rld"]])

  ## prepare the PCA
  mtx <- t(ret[["rld"]])
  vars <- apply(mtx, 2, var)
  sel  <- vars > 1e-26
  ret[["pca"]] <- prcomp(mtx[ , sel], scale.=TRUE)$x
  
  ret[["cntr_titles"]]        <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "ID")
  names(ret[["cntr_titles"]]) <- map_chr(ret[["config"]]$contrasts$contrast_list, `[[`, "title")
  ret[["cntr_titles"]]        <- ret[["cntr_titles"]][ ret[["cntr_titles"]] %in% names(ret[["cntr"]]) ]

  ret
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
