# Register server logic for the export-data module.
# Returns NULL; initializes file download handlers via bioshmods::fileExportServer.
.seapiper_server_export <- function(input, output, session, export_objects) {
  if(is.null(export_objects) || length(export_objects) == 0L) {
    return(NULL)
  }

  fileExportServer("fexp", export_objects)
}
