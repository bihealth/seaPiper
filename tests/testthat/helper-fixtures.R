# Build a small expression matrix fixture with stable gene/sample IDs.
make_fixture_exprs <- function() {
  exprs <- matrix(
    c(
      10, 12, 14, 16,
      20, 22, 24, 26,
      30, 29, 31, 33,
      40, 42, 41, 43
    ),
    nrow=4,
    byrow=TRUE
  )
  rownames(exprs) <- paste0("gene", seq_len(nrow(exprs)))
  colnames(exprs) <- paste0("sample", seq_len(ncol(exprs)))
  exprs
}

# Build annotation fixture and optionally include the PrimaryID column.
make_fixture_annot <- function(include_primary=TRUE) {
  annot <- data.frame(
    SYMBOL=paste0("G", seq_len(4)),
    stringsAsFactors=FALSE
  )
  rownames(annot) <- paste0("gene", seq_len(nrow(annot)))
  if(include_primary) {
    annot$PrimaryID <- rownames(annot)
  }
  annot
}

# Build one contrast table fixture; can remove list names and PrimaryID to test inference.
make_fixture_cntr <- function(named=TRUE, include_primary=FALSE) {
  cntr_tbl <- data.frame(
    log2FoldChange=c(1.2, -0.8, 0.4, -0.2),
    pvalue=c(0.001, 0.01, 0.4, 0.7),
    padj=c(0.01, 0.05, 0.5, 0.9),
    stringsAsFactors=FALSE
  )
  rownames(cntr_tbl) <- paste0("gene", seq_len(nrow(cntr_tbl)))

  if(include_primary) {
    cntr_tbl$PrimaryID <- rownames(cntr_tbl)
  }

  cntr <- list(cntr_a=cntr_tbl)
  if(!named) {
    names(cntr) <- NULL
  }
  cntr
}

# Build a simple covariate table keyed by sample ID.
make_fixture_covar <- function() {
  covar <- data.frame(
    ID=paste0("sample", seq_len(4)),
    group=c("A", "A", "B", "B"),
    stringsAsFactors=FALSE
  )
  rownames(covar) <- covar$ID
  covar
}

# Build a minimal but valid seaPiper config fixture.
make_fixture_config <- function(dataset_title="default") {
  list(
    organism=list(name="human", taxon="9606"),
    experiment=list(design_formula="~ group"),
    contrasts=list(
      contrast_list=list(
        list(ID="cntr_a", title="Contrast A")
      )
    ),
    tmod=list(
      databases=list(list(title="DB one")),
      sort_by=c("pvalue", "auc")
    ),
    filter=list(
      low_counts=10L,
      min_counts=5L,
      min_count_n=2L
    ),
    dataset_title=dataset_title
  )
}

# Assemble a full seaPiperData fixture using object-based constructor inputs.
make_fixture_spd <- function(
  include_primary_annot=TRUE,
  named_cntr=TRUE,
  include_primary_cntr=FALSE,
  title="default"
) {
  seapiperdata_from_objects(
    cntr=make_fixture_cntr(named=named_cntr, include_primary=include_primary_cntr),
    annot=make_fixture_annot(include_primary=include_primary_annot),
    exprs=make_fixture_exprs(),
    covar=make_fixture_covar(),
    config=make_fixture_config(dataset_title=title),
    primary_id="PrimaryID",
    title=title
  )
}

# Attach representative tmod structures to an existing seaPiperData fixture.
add_fixture_tmod <- function(spd) {
  spd$tmod_res <- list(
    default=list(
      cntr_a=list(
        DB1=list(
          pvalue=data.frame(score=0.1),
          auc=data.frame(score=0.2)
        )
      )
    )
  )
  spd$tmod_dbs <- list(
    default=list(
      DB1=list(dbobj=list(name="db1"))
    )
  )
  spd$tmod_map <- list(default=list(M1=c("gene1", "gene2")))
  spd$tmod_gl <- list(
    default=list(
      cntr_a=list(
        pvalue=setNames(c(1, 2), c("gene1", "gene2"))
      )
    )
  )
  spd
}

# Build a feature-flag list for mocked validate_seapiperdata() responses.
make_feature_flags <- function(value=FALSE) {
  list(
    info=value,
    gene_browser=value,
    heatmap=value,
    tmod=value,
    tmod_panel=value,
    disco=value,
    volcano=value,
    pca=value
  )
}
