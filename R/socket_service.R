
#' Disconnect socket-io server
#'
#' @return
#' @export
#'
init_socketio <- function() {
  io <- reticulate::import("socketio")
  e2$sio <- io$Client(logger = F, engineio_logger = F)
}


#' Disconnect socket-io server
#'
#' @return
#' @export
#'
disconnect_socketio <- function() {
  e2$sio$disconnect()
}

#' Send job progress to server
#'
#' @param message string
#' @export
#'
send_progress <- function(message = "message from R") {
  e2$sio$connect("http://127.0.0.1:9005")
  # e2$sio$connect("https://bmbls.bmi.osumc.edu/deepmaps/api/socket.io")
  # e2$sio$connect("https://bmbls.bmi.osumc.edu")
  e2$sio$emit("jobProgress", list(event = "jobProgress", data = message))
  # Don't know why the last emit will be ignored, could be async issue,
  # fixed this by emit a empty message
  Sys.sleep(0.1)
  e2$sio$emit("empty", "1")
  e2$sio$disconnect()
}

# i=1
# init_socketio()
# for(i in 1:10) {
#  Sys.sleep(0.5)
#  send_progress(i)
# }
