#' Echo working directory
#'
#' @param msg a character message
#'
#' @return a character message
#' @export
#'
print_wd <- function(msg = "Hello") {
  return(glue::glue("Message: {msg}
                    Working directory={getwd()}"))
}
