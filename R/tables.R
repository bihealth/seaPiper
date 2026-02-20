## create a table with the overview of the project
project_overview_table <- function(config, title) {
  
  tmp1 <- data.frame(
                 c("Organism", "Taxon ID", "Formula", "Contrasts", "Tmod dbs", 
                   "Low count filter (absolute)", "Low count filter (samples)"),
                 c(config$organism$name, 
                   config$organism$taxon,
                   config$experiment$design_formula,
                   length(config$contrasts$contrast_list),
                   paste(map_chr(config$tmod$databases, `[[`, "title"), collapse=", "),
                   sprintf("Removed genes with < %d total counts", config$filter$low_counts),
                   sprintf("Kept genes with at least %d counts in at least %d samples",
                           config$filter$min_counts, config$filter$min_count_n)
                   ))
  colnames(tmp1) <- c("Title:", title)
  tmp1
}

## create a table with contrast ids and titles
contrasts_overview_table <- function(config) {

  data.frame(ID=map_chr(config$contrasts$contrast_list, `[[`, "ID"),
             Title=map_chr(config$contrasts$contrast_list, `[[`, "title")
             )


}

## create a datatable with covariates
covariate_table <- function(covar) {

  tmp <- covar[ , !colnames(covar) %in% c("filename", "md5") ]

  datatable(tmp, extensions="Buttons",
              options=list(dom="Bfrtip", scrollX=TRUE, 
                           buttons=c("copy", "csv", "excel")))

}
