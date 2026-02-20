## server logic for the workflow info tab
.seapiper_server_info <- function(input, output, session, data, title) {
  observeEvent(input$select_pipeline, {
    ds <- input$select_pipeline
    if(!isTruthy(ds)) {
      return(NULL)
    }
                 
  output$project_overview   <- renderTable({ 
    project_overview_table(data[["config"]][[ds]], title) 
  })
  output$contrasts_overview <- renderTable({ 
    contrasts_overview_table(data[["config"]][[ds]]) })
  output$covariates         <- renderDT({ 
    covariate_table(data[["covar"]][[ds]]) 
    })
  output$session_info <- renderPrint(
                                     sessionInfo()
                                     )
  })
}
