test_that("merge_seapiperdata renames overlapping dataset IDs", {
  x <- make_fixture_spd(title="Dataset X")
  y <- make_fixture_spd(title="Dataset Y")

  merged <- merge_seapiperdata(x, y)

  expect_s3_class(merged, "seaPiperData")
  expect_setequal(names(merged$annot), c("default.x", "default.y"))
  expect_setequal(names(merged$cntr), c("default.x", "default.y"))
  expect_setequal(names(merged$config), c("default.x", "default.y"))
  expect_equal(merged$config$default.x$dataset_title, "Dataset X")
  expect_equal(merged$config$default.y$dataset_title, "Dataset Y")
})

test_that("merge_seapiperdata validates arguments", {
  spd <- make_fixture_spd()

  expect_error(
    merge_seapiperdata(list(), spd),
    "`x` must be a seaPiperData object"
  )

  expect_error(
    merge_seapiperdata(spd, list()),
    "`y` must be a seaPiperData object"
  )

  expect_error(
    merge_seapiperdata(spd, spd, suffixes=".x"),
    "must be a character vector of length 2"
  )

  expect_error(
    merge_seapiperdata(spd, spd, suffixes=c(".x", ".x")),
    "Duplicate dataset IDs"
  )
})
