#' Connect to socket-io server
#'
#' @param host string
#' @param port number
#' @return
#' @export
#'
connect_socketio <- function(host = "127.0.0.1", port = 9005) {
  io <- reticulate::import("socketio")
  #ws <- reticulate::import("websocket")
  #py_install("python-socketio[client]==4.6.1", pip=T)
  #py_install("python-socketio==4.6.1", pip=T)
  #py_install("websocket-client", pip=T)
  e2$sio <- io$Client(logger=F, engineio_logger=F)
  e2$sio$connect("http://127.0.0.1:9005")
  e2$sio$emit('jobProgress', list(event = "jobProgress", data = "Initializing computing service in R"))
}

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
  e2$sio$emit('jobProgress', list(event = "jobProgress", data = message))
}

#library(devtools)
