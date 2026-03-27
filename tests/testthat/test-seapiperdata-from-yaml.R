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
  expect_equal(spd$sample_id$default, "SampleID")
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

test_that("seapiperdata_from_yaml supports custom sample_id", {
  tmp_dir <- tempfile("sp_yaml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  covar <- make_fixture_covar()
  covar$sid <- covar$SampleID
  covar$SampleID <- NULL
  rownames(covar) <- covar$sid

  saveRDS(make_fixture_annot(), file.path(tmp_dir, "annot.rds"))
  saveRDS(make_fixture_cntr(), file.path(tmp_dir, "cntr.rds"))
  saveRDS(covar, file.path(tmp_dir, "covar.rds"))
  saveRDS(make_fixture_exprs(), file.path(tmp_dir, "rld.rds"))

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    sample_id: sid",
      "    annot: annot.rds",
      "    cntr: cntr.rds",
      "    covar: covar.rds",
      "    rld: rld.rds"
    ),
    file.path(tmp_dir, "custom_sample_id.yaml")
  )

  spd <- seapiperdata_from_yaml(yaml_file)
  expect_equal(spd$sample_id$default, "sid")
  expect_equal(rownames(spd$covar$default), covar$sid)
})

test_that("seapiperdata_from_yaml validates user-provided config explicitly", {
  tmp_dir <- tempfile("sp_yaml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  saveRDS(make_fixture_annot(), file.path(tmp_dir, "annot.rds"))
  saveRDS(make_fixture_cntr(), file.path(tmp_dir, "cntr.rds"))
  saveRDS(make_fixture_covar(), file.path(tmp_dir, "covar.rds"))
  saveRDS(make_fixture_exprs(), file.path(tmp_dir, "rld.rds"))

  bad_config <- list(dataset_title="broken")
  saveRDS(bad_config, file.path(tmp_dir, "config_missing_fields.rds"))

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    annot: annot.rds",
      "    cntr: cntr.rds",
      "    covar: covar.rds",
      "    rld: rld.rds",
      "    config: config_missing_fields.rds"
    ),
    file.path(tmp_dir, "bad_config.yaml")
  )

  expect_error(
    seapiperdata_from_yaml(yaml_file),
    "config for dataset `default` is missing required field\\(s\\): organism, experiment, contrasts, tmod, filter"
  )

  bad_config <- make_fixture_config()
  bad_config$filter$min_count_n <- NULL
  saveRDS(bad_config, file.path(tmp_dir, "config_missing_filter_field.rds"))

  yaml_file <- write_yaml_fixture(
    c(
      "datasets:",
      "  default:",
      "    annot: annot.rds",
      "    cntr: cntr.rds",
      "    covar: covar.rds",
      "    rld: rld.rds",
      "    config: config_missing_filter_field.rds"
    ),
    file.path(tmp_dir, "bad_filter_config.yaml")
  )

  expect_error(
    seapiperdata_from_yaml(yaml_file),
    "config for dataset `default`: `filter` is missing required field\\(s\\): min_count_n"
  )
})
