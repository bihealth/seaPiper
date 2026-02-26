# Write YAML fixture content to disk and return the written file path.
write_yaml_fixture <- function(lines, path) {
  writeLines(lines, con=path)
  path
}

test_that("seapiperdata_from_yaml loads relative RDS paths and computes pca", {
  tmp_dir <- tempfile("sp_yaml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  saveRDS(make_fixture_annot(), file.path(tmp_dir, "annot.rds"))
  saveRDS(make_fixture_cntr(), file.path(tmp_dir, "cntr.rds"))
  saveRDS(make_fixture_covar(), file.path(tmp_dir, "covar.rds"))
  saveRDS(make_fixture_exprs(), file.path(tmp_dir, "rld.rds"))
  saveRDS(make_fixture_config(dataset_title="YAML dataset"), file.path(tmp_dir, "config.rds"))

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    annot: annot.rds",
      "    cntr: cntr.rds",
      "    covar: covar.rds",
      "    rld: rld.rds",
      "    config: config.rds"
    ),
    file.path(tmp_dir, "seapiper_data.yaml")
  )

  spd <- seapiperdata_from_yaml(yaml_file)

  expect_s3_class(spd, "seaPiperData")
  expect_true(is.matrix(spd$pca$default))
  expect_equal(spd$config$default$dataset_title, "YAML dataset")
})

test_that("seapiperdata_from_yaml validates YAML structure and required keys", {
  tmp_dir <- tempfile("sp_yaml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  yaml_file <- write_yaml_fixture(
    c("not_datasets:", "  default:", "    x: 1"),
    file.path(tmp_dir, "invalid.yaml")
  )

  expect_error(
    seapiperdata_from_yaml(yaml_file),
    "top-level `datasets` mapping"
  )

  saveRDS(make_fixture_annot(), file.path(tmp_dir, "annot.rds"))
  saveRDS(make_fixture_cntr(), file.path(tmp_dir, "cntr.rds"))
  saveRDS(make_fixture_exprs(), file.path(tmp_dir, "rld.rds"))

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    annot: annot.rds",
      "    cntr: cntr.rds",
      "    rld: rld.rds"
    ),
    file.path(tmp_dir, "missing_covar.yaml")
  )

  expect_error(
    seapiperdata_from_yaml(yaml_file),
    "`covar` is required"
  )
})

test_that("seapiperdata_from_yaml errors when referenced file is missing", {
  tmp_dir <- tempfile("sp_yaml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    covar: covar.rds"
    ),
    file.path(tmp_dir, "missing_file.yaml")
  )

  expect_error(
    seapiperdata_from_yaml(yaml_file),
    "file for `covar` not found"
  )
})
