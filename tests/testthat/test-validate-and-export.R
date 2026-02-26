test_that("validate_seapiperdata sets feature flags for minimal and full inputs", {
  spd_min <- make_fixture_spd()
  res_min <- seaPiper:::validate_seapiperdata(spd_min)

  expect_true(res_min$features$gene_browser)
  expect_true(res_min$features$volcano)
  expect_true(res_min$features$disco)
  expect_true(res_min$features$pca)
  expect_true(res_min$features$info)
  expect_false(res_min$features$tmod)
  expect_false(res_min$features$tmod_panel)
  expect_setequal(names(res_min$missing), c("tmod", "tmod_panel"))

  spd_full <- add_fixture_tmod(spd_min)
  res_full <- seaPiper:::validate_seapiperdata(spd_full)

  expect_true(res_full$features$tmod)
  expect_true(res_full$features$tmod_panel)
  expect_false("tmod" %in% names(res_full$missing))
  expect_false("tmod_panel" %in% names(res_full$missing))
})

test_that("export dataframe payload helpers normalize supported shapes", {
  mat <- matrix(1:4, nrow=2)
  payload <- seaPiper:::.as_export_dataframe_payload(mat)
  expect_s3_class(payload, "data.frame")

  lvl1 <- list(a=data.frame(x=1), b=matrix(2:3, nrow=1))
  payload <- seaPiper:::.as_export_dataframe_payload(lvl1)
  expect_true(is.list(payload))
  expect_true(all(vapply(payload, is.data.frame, logical(1))))

  lvl2 <- list(
    ds1=list(a=data.frame(x=1), b=data.frame(x=2)),
    ds2=list(a=matrix(3:4, nrow=1))
  )
  payload <- seaPiper:::.as_export_dataframe_payload(lvl2)
  expect_true(is.list(payload))
  expect_true(all(vapply(payload, is.list, logical(1))))
  expect_true(all(vapply(unlist(payload, recursive=FALSE), is.data.frame, logical(1))))
})

test_that("collapse nested list of dataframes behaves and validates input", {
  nested <- list(
    DB1=list(
      pvalue=data.frame(score=1),
      auc=data.frame(score=2)
    ),
    DB2=list(
      pvalue=data.frame(score=3)
    )
  )

  collapsed <- seaPiper:::.collapse_nested_ldf(nested, c("Database", "Sorting"))
  expect_s3_class(collapsed, "data.frame")
  expect_equal(colnames(collapsed), c("Database", "Sorting", "score"))
  expect_equal(nrow(collapsed), 3L)
  expect_setequal(unique(collapsed$Database), c("DB1", "DB2"))

  expect_error(
    seaPiper:::.collapse_nested_ldf(list(), c("Database")),
    "non-empty list"
  )
  expect_error(
    seaPiper:::.collapse_nested_ldf(nested, character(0)),
    "No more column names"
  )
})

test_that("tmod result export conversion returns dataframe tables", {
  tmod_res <- list(
    default=list(
      cntr_a=list(
        DB1=list(
          pvalue=data.frame(score=1),
          auc=data.frame(score=2)
        )
      )
    )
  )

  capture.output(value <- seaPiper:::.tmod_res_as_lldf(tmod_res))
  expect_s3_class(value$default$cntr_a, "data.frame")
  expect_equal(colnames(value$default$cntr_a), c("Database", "Sorting", "score"))
})

test_that("build_export_objects includes expected core exports", {
  spd <- make_fixture_spd()
  specs <- seaPiper:::.build_export_objects(spd)

  titles <- vapply(specs, `[[`, character(1), "title")
  expect_true("Differential Expression Contrasts" %in% titles)
  expect_true("Gene Annotation" %in% titles)
  expect_true("Covariates" %in% titles)
  expect_true("Expression Matrix (rlog)" %in% titles)
  expect_true("PCA Coordinates" %in% titles)
  expect_true("Pipeline Configuration (Raw Object)" %in% titles)
  expect_true("Complete seaPiperData Object" %in% titles)
  expect_false("Tmod Results" %in% titles)
})
