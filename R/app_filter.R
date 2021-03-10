#' Simple request logger
#'
#' @param req plumber request environment
#' @param res plumber response environment
#'
#' @export
log_request <- function(req, res) {
  cat(
    as.character(Sys.time()), "-",
    req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  plumber::forward()
}
