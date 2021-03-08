#' Print working directory and a message
#'
#' @param msg a character message
#'
#' @return a character message
#'
print_wd <- function(msg = "Hello") {
  return(glue::glue("Message: {msg},\\
                    Working directory={getwd()}"))
}
