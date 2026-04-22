# Global logging switch for package-internal tracing.
# By default, logging remains enabled to preserve current user-facing messages.
# Disable with: options(seapiper.debug = FALSE)
.seapiper_log_enabled <- function() {
  opt <- getOption("seapiper.debug", NULL)
  if(is.null(opt)) {
    return(TRUE)
  }

  isTRUE(opt)
}

# Internal logging helper with optional file/function prefix.
.seapiper_log <- function(..., .prefix=NULL) {
  if(!.seapiper_log_enabled()) {
    return(invisible(NULL))
  }

  msg <- paste(..., collapse="")
  if(!is.null(.prefix) && nzchar(.prefix)) {
    msg <- paste0("[", .prefix, "] ", msg)
  }

  message(msg)
  invisible(NULL)
}
