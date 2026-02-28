## server logic for gene browser modules
.seapiper_server_gene <- function(input, output, session, data, gene_id, palettes) {
  geneBrowserTableServer(id="geneT",
    cntr=data[["cntr"]],
    annot=data[["annot"]],
    annot_linkout=data[["annot_linkout"]],
    gene_id=gene_id)

  observeEvent(gene_id$id, {
    updateTabItems(session, "navid", "gene_browser")
  })

  geneBrowserPlotServer(id="geneP",
                        gene_id=gene_id,
                        covar=data[["covar"]],
                        exprs=data[["rld"]], annot=data[["annot"]], 
                        annot_linkout=data[["annot_linkout"]],
                        cntr=data[["cntr"]],
                        exprs_label = "Regularized log transformed expression (rlog)",
                        palettes=palettes
  )
}
