test_that("seapiper enforces bioshmods version requirement", {
  local_mocked_bindings(
    # Mock installed bioshmods version as older than required.
    packageVersion=function(...) numeric_version("0.0.0"),
    # Force version comparison to report installed < required.
    compareVersion=function(...) -1L,
    .package="utils"
  )

  data <- structure(list(), class="seaPiperData")
  expect_error(
    seapiper(data),
    "requires bioshmods >="
  )
})

test_that("seapiper validates input class", {
  local_mocked_bindings(
    # Mock installed bioshmods version as sufficiently new.
    packageVersion=function(...) numeric_version("9.9.9"),
    # Force version comparison to pass the version gate.
    compareVersion=function(...) 1L,
    .package="utils"
  )

  expect_error(
    seapiper(list()),
    "`data` must be a seaPiperData object"
  )
})

test_that("seapiper can return a shiny app object with mocked internals", {
  local_mocked_bindings(
    # Mock installed bioshmods version as sufficiently new.
    packageVersion=function(...) numeric_version("9.9.9"),
    # Force version comparison to pass the version gate.
    compareVersion=function(...) 1L,
    .package="utils"
  )

  local_mocked_bindings(
    # Return a minimal validation response with all feature flags disabled.
    validate_seapiperdata=function(...) {
      list(
        features=make_feature_flags(FALSE),
        missing=list(),
        present=list()
      )
    },
    # Return no export objects so export tab stays disabled.
    .build_export_objects=function(...) list(),
    # Build a valid dashboard header for shiny.appobj construction.
    .pipeline_dashboard_header=function(title) shinydashboard::dashboardHeader(title=title),
    # Build an empty but valid dashboard sidebar.
    .pipeline_dashboard_sidebar=function(...) shinydashboard::dashboardSidebar(),
    # Build an empty but valid dashboard body.
    .pipeline_dashboard_body=function(...) shinydashboard::dashboardBody(),
    .package="seaPiper"
  )

  data <- structure(list(), class="seaPiperData")
  app <- seapiper(data, title="App title")
  expect_s3_class(app, "shiny.appobj")
})

test_that("seapiper reports disabled modules when inputs are missing", {
  local_mocked_bindings(
    # Mock installed bioshmods version as sufficiently new.
    packageVersion=function(...) numeric_version("9.9.9"),
    # Force version comparison to pass the version gate.
    compareVersion=function(...) 1L,
    .package="utils"
  )

  local_mocked_bindings(
    # Return validation with one missing module dependency to trigger message().
    validate_seapiperdata=function(...) {
      list(
        features=make_feature_flags(FALSE),
        missing=list(volcano="cntr"),
        present=list()
      )
    },
    # Return no export objects so export tab stays disabled.
    .build_export_objects=function(...) list(),
    # Build a valid dashboard header for shiny.appobj construction.
    .pipeline_dashboard_header=function(...) shinydashboard::dashboardHeader(),
    # Build an empty but valid dashboard sidebar.
    .pipeline_dashboard_sidebar=function(...) shinydashboard::dashboardSidebar(),
    # Build an empty but valid dashboard body.
    .pipeline_dashboard_body=function(...) shinydashboard::dashboardBody(),
    .package="seaPiper"
  )

  data <- structure(list(), class="seaPiperData")
  expect_message(
    seapiper(data),
    "disabled modules"
  )
})
