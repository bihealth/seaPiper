test_that("merge_seapiperdata supports explicit dataset names", {
  x <- make_fixture_spd(title="Dataset X")
  y <- make_fixture_spd(title="Dataset Y")

  merged <- merge_seapiperdata(x, y, names=c("one", "two"))

  expect_s3_class(merged, "seaPiperData")
  expect_setequal(names(merged$annot), c("one", "two"))
  expect_setequal(names(merged$cntr), c("one", "two"))
  expect_setequal(names(merged$config), c("one", "two"))
  expect_equal(merged$config$one$dataset_title, "Dataset X")
  expect_equal(merged$config$two$dataset_title, "Dataset Y")
})

test_that("merge_seapiperdata supports more than two inputs", {
  x <- make_fixture_spd(title="Dataset X")
  y <- make_fixture_spd(title="Dataset Y")
  z <- make_fixture_spd(title="Dataset Z")

  merged <- merge_seapiperdata(x, y, z, names=c("one", "two", "three"))

  expect_s3_class(merged, "seaPiperData")
  expect_setequal(names(merged$annot), c("one", "two", "three"))
  expect_setequal(names(merged$cntr), c("one", "two", "three"))
  expect_setequal(names(merged$config), c("one", "two", "three"))
  expect_equal(merged$config$three$dataset_title, "Dataset Z")
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
    merge_seapiperdata(spd, spd, names="x"),
    "must be a character vector of length 2"
  )

  expect_error(
    merge_seapiperdata(spd, spd),
    "Duplicate dataset IDs"
  )
})

test_that("name<- renames dataset IDs across all fields", {
  x <- make_fixture_spd(title="Dataset X")
  y <- make_fixture_spd(title="Dataset Y")
  merged <- merge_seapiperdata(x, y, names=c("x", "y"))

  name(merged) <- c("one", "two")

  expect_setequal(names(merged$annot), c("one", "two"))
  expect_setequal(names(merged$cntr), c("one", "two"))
  expect_setequal(names(merged$config), c("one", "two"))
})
