## server logic for gene browser modules
.seapiper_server_gene <- function(input, output, session, data, gene_id) {
  geneBrowserTableServer("geneT", data[["cntr"]], data[["annot"]], 
    annot_linkout=data[["annot_linkout"]],
    gene_id=gene_id)

  observeEvent(gene_id$id, {
    updateTabItems(session, "navid", "gene_browser")
  })

  geneBrowserPlotServer("geneP", gene_id, covar=data[["covar"]], 
                        exprs=data[["rld"]], annot=data[["annot"]], 
                        annot_linkout=data[["annot_linkout"]],
                        cntr=data[["cntr"]],
                        exprs_label = "Regularized log transformed expression (rlog)"
  )
}
