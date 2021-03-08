#' Start IRIS3 R api server
#' @import plumber
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @return NULL
#' @examples
#' \dontrun{
#' start_server(port = 8000)
#' }
#' @export
start_server <-
  function(port = 8000,
           host = "127.0.0.1") {
    message(paste("Starting IRIS3 R server at host =", host, ", port =", port))
    tryCatch(
      {
        plumber::plumb(dir = as.character(system.file("endpoints",
          package = "iris3api"
        ))) %>%
          plumber::pr_set_docs(FALSE) %>% # Disable swagger
          plumber::pr_run(
            host = host,
            port = as.numeric(port)
          )
      },
      interrupt = function(i) {
        message("Server Exited")
      }
    )
    return(NULL)
  }
