#' Start IRIS3 R api server
#' @import plumber
#' @param port default: 8000
#' @param host default: "127.0.0.1"
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
          plumber::pr_hook("preroute", function(req, res) {
            tictoc::tic()
            logger::log_info(
              glue(
               '{convert_empty(req$REMOTE_ADDR)} \\
              "{convert_empty(req$HTTP_USER_AGENT)}"\\
               {convert_empty(req$REQUEST_METHOD)}\\
               {convert_empty(req$PATH_INFO)}\\
               {convert_empty(res$status)}'
              )
            )
          }) %>%
          plumber::pr_hook("postroute", function(req, res) {
            end <- tictoc::toc(quiet = TRUE)
            logger::log_info(
              glue(
                '
               {convert_empty(req$REQUEST_METHOD)} \\
               {convert_empty(req$PATH_INFO)} \\
               {convert_empty(res$status)} \\
               Time elapsed: \\
               {round(end$toc-end$tic,digits=5)}s'
              )
            )
          }) %>%
          plumber::pr_set_docs(FALSE) %>%
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
