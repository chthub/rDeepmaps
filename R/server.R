#' Start IRIS3 R api server
#' @importFrom jsonlite fromJSON
#' @import plumber
#' @importFrom methods new
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @return object
#' @export
start_server <-
  function(port = NULL,
           host = NULL) {
    if (is.null(port)) {
      port <- sample(8000:9999, 1)
    }
    if (is.null(host)) {
      host <- "127.0.0.1"
    }

    # Define the API
    app <- plumber::pr()
    controller <- plumber::pr("R/controller.R")

    message(paste("Starting IRIS3 R server at host =", host, ", port =", port))
    tryCatch(
      {
        app %>%
          pr_set_docs(FALSE) %>% # Disable swagger
          plumber::pr_mount("/", controller) %>%
          plumber::pr_run(
            host = host,
            port = port
          )
      },
      interrupt = function(i) {
        message("Server Exited")
      }
    )
    return(app)
  }


#' DEV tool: auto style and lint package
#'
#' @return NULL
#' @export
#'
style_and_lint <- function() {
  styler::style_pkg()
  lintr::lint_package()
}
