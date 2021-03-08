#' Echo the parameter that was sent in
#' @name Echo
#' @param msg The message to echo back.
#'
#* @get /echo
function(msg = "") {
  return(print_wd(msg))
}
