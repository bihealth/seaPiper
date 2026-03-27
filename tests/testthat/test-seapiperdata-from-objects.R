test_that("seapiperdata_from_objects builds expected seaPiperData shape", {
  spd <- make_fixture_spd()

  expect_s3_class(spd, "seaPiperData")
  expect_true("default" %in% names(spd$annot))
  expect_true("default" %in% names(spd$cntr))
  expect_true("default" %in% names(spd$covar))
  expect_true("default" %in% names(spd$rld))
  expect_true("default" %in% names(spd$pca))
  expect_true("default" %in% names(spd$config))

  expect_true("PrimaryID" %in% colnames(spd$annot$default))
  expect_equal(spd$config$default$dataset_title, "default")
  expect_true(is.matrix(spd$rld$default))
  expect_true(is.matrix(spd$pca$default))
  expect_equal(spd$sample_id$default, "SampleID")
})

test_that("seapiperdata_from_objects infers primary IDs when missing", {
  spd <- make_fixture_spd(
    include_primary_annot=FALSE,
    include_primary_cntr=FALSE
  )

  expect_true("PrimaryID" %in% colnames(spd$annot$default))
  expect_true("PrimaryID" %in% colnames(spd$cntr$default$cntr_a))
  expect_identical(spd$annot$default$PrimaryID, rownames(spd$annot$default))
})

test_that("seapiperdata_from_objects assigns contrast names when unnamed", {
  spd <- seapiperdata_from_objects(
    cntr=make_fixture_cntr(named=FALSE),
    annot=make_fixture_annot(),
    exprs=make_fixture_exprs(),
    covar=make_fixture_covar(),
    config=NULL
  )

  expect_equal(names(spd$cntr$default), "contrast_1")
  expect_equal(unname(spd$cntr_titles$default), "contrast_1")
})

test_that("seapiperdata_from_objects validates expression matrix names", {
  exprs <- make_fixture_exprs()
  colnames(exprs) <- NULL

  expect_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=exprs
    ),
    "must have column names"
  )

  exprs <- make_fixture_exprs()
  rownames(exprs) <- NULL

  expect_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=exprs
    ),
    "must have row names"
  )
})

test_that("seapiperdata_from_objects honors custom sample_id and validates covariates", {
  covar <- make_fixture_covar()
  covar$sample_label <- covar$SampleID
  covar$SampleID <- NULL
  rownames(covar) <- covar$sample_label

  spd <- seapiperdata_from_objects(
    cntr=make_fixture_cntr(),
    annot=make_fixture_annot(),
    exprs=make_fixture_exprs(),
    covar=covar,
    sample_id="sample_label"
  )

  expect_equal(spd$sample_id$default, "sample_label")
  expect_equal(rownames(spd$covar$default), covar$sample_label)

  expect_no_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=make_fixture_exprs(),
      covar=make_fixture_covar()
    )
  )

  bad_covar <- make_fixture_covar()
  bad_covar$SampleID <- NULL
  expect_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=make_fixture_exprs(),
      covar=bad_covar
    ),
    "must contain sample ID column"
  )
})

test_that("seapiperdata_from_objects validates user-provided config explicitly", {
  expect_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=make_fixture_exprs(),
      covar=make_fixture_covar(),
      config=list(dataset_title="broken")
    ),
    "missing required field\\(s\\): organism, experiment, contrasts, tmod, filter"
  )

  bad_config <- make_fixture_config()
  bad_config$contrasts$contrast_list <- list(list(ID="cntr_a"))

  expect_error(
    seapiperdata_from_objects(
      cntr=make_fixture_cntr(),
      annot=make_fixture_annot(),
      exprs=make_fixture_exprs(),
      covar=make_fixture_covar(),
      config=bad_config
    ),
    "each entry in `contrasts\\$contrast_list` must contain non-empty `ID` and `title` strings"
  )
})
