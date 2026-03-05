## check whether a list field contains usable data entries
.has_list_data <- function(x, strict=FALSE) {
  if(!is.list(x) || length(x) == 0L) {
    return(FALSE)
  }

  has_data <- vapply(x, function(.x) !is.null(.x) && length(.x) > 0L, logical(1))
  if(isTRUE(strict)) {
    all(has_data)
  } else {
    any(has_data)
  }
}
