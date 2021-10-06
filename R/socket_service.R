
#' Disconnect socket-io server
#'
#' @return
#' @export
#'
init_socketio <- function() {
  #message(paste("Starting socket.io"))
  #io <- reticulate::import("socketio")
  #e2$sio <- io$Client(logger = F, engineio_logger = F)
  #e2$sio$connect("http://127.0.0.1:9005")
  #e2$sio$connect("https://bmblx.bmi.osumc.edu")
  #e2$sio$connect("http://10.82.14.183:9005")
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
send_progress <- function(message = "message from R1") {
  #e2$sio$emit('jobProgress', list(event = "jobProgress", data = message))
  # Don't know why the last emit will be ignored, could be async issue,
  # fixed this by emit a empty message
  #Sys.sleep(0.1)
  #e2$sio$emit('empty', "1")
}

#init_socketio()
#
# i=1
# for(i in 1:10) {
#  Sys.sleep(0.2)
#  send_progress(i)
# }

