## server logic for tmod-related modules
.seapiper_server_tmod <- function(input, output, session, data, gs_id, gene_id, enable_panel=FALSE) {
  tmodBrowserPlotServer("tmodP", gs_id, 
                                    tmod_dbs=data[["tmod_dbs"]], 
                                    cntr    =data[["cntr"]], 
                                    tmod_map=data[["tmod_map"]], 
                                    tmod_gl =data[["tmod_gl"]], 
                                    annot   =data[["annot"]],
                                    tmod_res=data[["tmod_res"]],
                                    gene_id=gene_id)

  tmodBrowserTableServer("tmodT", data[["tmod_res"]], gs_id=gs_id, 
                         multilevel=TRUE, tmod_dbs=data[["tmod_dbs"]])

  observeEvent(gs_id$id, { 
    updateTabItems(session, "navid", "tmod_browser")
  })

  if(isTRUE(enable_panel)) {
    tmodPanelPlotServer("panelP", cntr    =data[["cntr"]], 
                                  tmod_res=data[["tmod_res"]],
                                  tmod_dbs=data[["tmod_dbs"]], 
                                  tmod_map=data[["tmod_map"]], 
                                  annot   =data[["annot"]],
                                  gs_id=gs_id)
  }
}
