test_that("project_overview_table returns expected summary rows", {
  config <- make_fixture_config(dataset_title="Table dataset")
  tbl <- project_overview_table(config, title="Pipeline title")

  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 7L)
  expect_equal(colnames(tbl), c("Title:", "Pipeline title"))
  expect_true(any(grepl("Organism", tbl[[1]])))
})

test_that("contrasts_overview_table returns configured contrast IDs and titles", {
  config <- make_fixture_config()
  tbl <- contrasts_overview_table(config)

  expect_s3_class(tbl, "data.frame")
  expect_equal(colnames(tbl), c("ID", "Title"))
  expect_equal(tbl$ID, "cntr_a")
  expect_equal(tbl$Title, "Contrast A")
})

test_that("covariate_table excludes filename and md5 fields", {
  covar <- make_fixture_covar()
  covar$filename <- "fileA"
  covar$md5 <- "abc123"

  widget <- covariate_table(covar)
  expect_s3_class(widget, "datatables")
  expect_false("filename" %in% colnames(widget$x$data))
  expect_false("md5" %in% colnames(widget$x$data))
})
