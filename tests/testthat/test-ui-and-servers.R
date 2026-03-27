test_that("get_dataset_ids returns first available named list section", {
  data <- list(
    covar=list(ds_cov=data.frame(ID="s1", row.names="s1")),
    config=list(ds_cfg=list(dataset_title="cfg"))
  )

  ids <- seaPiper:::.get_dataset_ids(data)
  expect_identical(ids, "ds_cov")

  ids <- seaPiper:::.get_dataset_ids(list(a=1))
  expect_equal(ids, character(0))
})

test_that("expression modules use a single explicit sample-id column", {
  expect_equal(
    seaPiper:::.resolve_sample_id_col(list(default="SampleID")),
    "SampleID"
  )

  expect_error(
    seaPiper:::.resolve_sample_id_col(list(ds1="SampleID", ds2="sid")),
    "same `sample_id`"
  )

  expect_no_error(
    seaPiper:::.normalize_expression_covar(
      covar=list(default=make_fixture_covar()),
      sample_id_col="SampleID"
    )
  )

  bad_covar <- make_fixture_covar()
  bad_covar$SampleID <- NULL
  expect_error(
    seaPiper:::.normalize_expression_covar(
      covar=list(default=bad_covar),
      sample_id_col="SampleID"
    ),
    "must contain sample ID column"
  )
})

test_that("sidebar always includes Help and conditionally includes module tabs", {
  features <- list(
    gene_browser=TRUE,
    volcano=FALSE,
    disco=FALSE,
    tmod=FALSE,
    tmod_panel=FALSE,
    pca=FALSE,
    file_export=FALSE,
    info=FALSE
  )

  sidebar <- seaPiper:::.pipeline_dashboard_sidebar(features)
  html <- as.character(sidebar)

  expect_match(html, "Help")
  expect_match(html, "Gene browser")
  expect_no_match(html, "Volcano plots")
})

test_that("pipeline body can be built with minimal features", {
  data <- list(
    covar=list(default=make_fixture_covar()),
    config=list(default=make_fixture_config(dataset_title="Dataset One")),
    cntr=NULL,
    annot=NULL,
    rld=NULL,
    pca=NULL,
    tmod_res=NULL,
    tmod_dbs=NULL,
    tmod_map=NULL,
    tmod_gl=NULL,
    cntr_titles=NULL
  )

  features <- list(
    gene_browser=FALSE,
    volcano=FALSE,
    disco=FALSE,
    tmod=FALSE,
    tmod_panel=FALSE,
    pca=FALSE,
    file_export=FALSE,
    info=FALSE
  )

  body <- seaPiper:::.pipeline_dashboard_body(data, "Test title", features, export_objects=list())
  html <- as.character(body)
  expect_match(html, "help")
})

test_that("server misc wrapper dispatches only enabled modules", {
  calls <- new.env(parent=emptyenv())
  calls$disco <- 0L
  calls$volcano <- 0L
  calls$heatmap <- 0L
  calls$pca <- 0L
  calls$selected_ids_volcano <- NULL
  calls$selected_ids_heatmap <- NULL
  calls$sample_id_col_heatmap <- NULL

  selected_ids <- list(ids=c("gene1", "gene2"))

  local_mocked_bindings(
    # Count heatmap module dispatches and capture shared selected_ids wiring.
    heatmapServer=function(...) {
      calls$heatmap <- calls$heatmap + 1L
      args <- list(...)
      calls$selected_ids_heatmap <- args$selected_ids
      calls$sample_id_col_heatmap <- args$sample_id_col
    },
    .package="bioshmods"
  )

  local_mocked_bindings(
    # Count disco module dispatches.
    discoServer=function(...) calls$disco <- calls$disco + 1L,
    # Count volcano module dispatches and capture shared selected_ids wiring.
    volcanoServer=function(...) {
      calls$volcano <- calls$volcano + 1L
      args <- list(...)
      calls$selected_ids_volcano <- args$selected_ids
    },
    # Count PCA module dispatches.
    pcaServer=function(...) calls$pca <- calls$pca + 1L,
    .package="seaPiper"
  )

  data <- list(cntr=list(), annot=list(), pca=list(), covar=list(), rld=list())
  seaPiper:::.seapiper_server_misc(
    input=NULL,
    output=NULL,
    session=NULL,
    data=data,
    gene_id=list(),
    selected_ids=selected_ids,
    enable_disco=TRUE,
    enable_volcano=TRUE,
    enable_heatmap=TRUE,
    enable_pca=TRUE,
    sample_id_col="SampleID"
  )

  expect_equal(calls$disco, 1L)
  expect_equal(calls$volcano, 1L)
  expect_equal(calls$heatmap, 1L)
  expect_equal(calls$pca, 1L)
  expect_identical(calls$selected_ids_volcano, selected_ids)
  expect_identical(calls$selected_ids_heatmap, selected_ids)
  expect_identical(calls$sample_id_col_heatmap, "SampleID")
})

test_that("server export wrapper only registers when payload is present", {
  calls <- new.env(parent=emptyenv())
  calls$n <- 0L

  local_mocked_bindings(
    # Count export server registrations.
    fileExportServer=function(...) calls$n <- calls$n + 1L,
    .package="seaPiper"
  )

  seaPiper:::.seapiper_server_export(NULL, NULL, NULL, NULL)
  seaPiper:::.seapiper_server_export(NULL, NULL, NULL, list())
  expect_equal(calls$n, 0L)

  export_objects <- list(
    list(
      data=data.frame(x=1),
      title="X",
      description="desc",
      data_type="dataframes"
    )
  )
  seaPiper:::.seapiper_server_export(NULL, NULL, NULL, export_objects)
  expect_equal(calls$n, 1L)
})
