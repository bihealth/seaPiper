## server logic for disco, volcano, and pca modules
.seapiper_server_misc <- function(input, output, session, data, gene_id, enable_disco=FALSE,
                                  enable_volcano=FALSE, enable_pca=FALSE) {
  if(isTRUE(enable_disco)) {
    discoServer("disco", data[["cntr"]], data[["annot"]], gene_id=gene_id)
  }

  if(isTRUE(enable_volcano)) {
    volcanoServer("volcano", data[["cntr"]], annot=data[["annot"]], gene_id=gene_id)
  }

  if(isTRUE(enable_pca)) {
    pcaServer("pca", data[["pca"]], data[["covar"]])
  }
}
