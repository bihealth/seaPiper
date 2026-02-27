## validate a seaPiperData object and report missing inputs
validate_seapiperdata <- function(data) {
  has_list_data <- function(x) {
    is.list(x) && length(x) > 0 &&
      any(vapply(x, function(.x) !is.null(.x) && length(.x) > 0, logical(1)))
  }

  present <- list(
    cntr      = has_list_data(data[["cntr"]]),
    annot     = has_list_data(data[["annot"]]),
    covar     = has_list_data(data[["covar"]]),
    rld       = has_list_data(data[["rld"]]),
    pca       = has_list_data(data[["pca"]]),
    config    = has_list_data(data[["config"]]),
    tmod_res  = has_list_data(data[["tmod_res"]]),
    tmod_dbs  = has_list_data(data[["tmod_dbs"]]),
    tmod_map  = has_list_data(data[["tmod_map"]]),
    tmod_gl   = has_list_data(data[["tmod_gl"]])
  )

  features <- list(
    covar        = present$covar,
    gene_browser = present$cntr && present$annot && present$covar && present$rld,
    volcano      = present$cntr,
    heatmap      = present$annot && present$rld,
    disco        = present$cntr && present$annot,
    pca          = present$pca && present$covar,
    info         = present$config && present$covar,
    tmod         = present$tmod_res && present$tmod_dbs && present$tmod_map && present$tmod_gl,
    tmod_panel   = present$tmod_res && present$tmod_dbs && present$tmod_map && present$annot
  )

  missing <- list(
    gene_browser = c("cntr", "annot", "covar", "rld")[!c(present$cntr, present$annot, present$covar, present$rld)],
    volcano      = c("cntr")[!present$cntr],
    heatmap      = c("annot", "rld")[!c(present$annot, present$rld)],
    disco        = c("cntr", "annot")[!c(present$cntr, present$annot)],
    pca          = c("pca", "covar")[!c(present$pca, present$covar)],
    info         = c("config", "covar")[!c(present$config, present$covar)],
    tmod         = c("tmod_res", "tmod_dbs", "tmod_map", "tmod_gl")[
                    !c(present$tmod_res, present$tmod_dbs, present$tmod_map, present$tmod_gl)],
    tmod_panel   = c("tmod_res", "tmod_dbs", "tmod_map", "annot")[
                    !c(present$tmod_res, present$tmod_dbs, present$tmod_map, present$annot)]
  )

  missing <- missing[vapply(missing, length, integer(1)) > 0]

  list(features=features, missing=missing, present=present)
}
