
#' Disconnect socket-io server
#'
#' @param host string
#' @param port number
#' @return
#' @export
#'
disconnect_socketio <- function() {
  e2$sio$disconnect()
}

#' Send job progress to server
#'
#' @param message
#' @return
#' @export
#'
send_progress <- function(message = "message from R") {
  io <- reticulate::import("socketio")
  e2$sio <- io$Client(logger = F, engineio_logger = F)
  e2$sio$connect("http://127.0.0.1:9005")
  e2$sio$emit('jobProgress', list(event = "jobProgress", data = message))
  # Don't know why the last emit will be ignored, could be async issue,
  # fixed this by emit a empty message
  Sys.sleep(0.1)
  e2$sio$emit('empty', "1")
  #e2$sio$disconnect()
}

#library(reticulate)
#py_install("python-socketio[client]==4.6.1", pip=T)
#py_install("python-socketio==4.6.1", pip=T)
#py_install("websocket-client", pip=T)
