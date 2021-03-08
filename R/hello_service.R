#' Print working directory and a message
#'
#' @param msg a character message
#' @return a character message
#' @export
#' @examples
#' print_wd(msg = "Hi")
print_wd <- function(msg = "Hello") {
  return(glue::glue("Message: {msg},\\
                    Working directory={getwd()}"))
}

#' Echo back a message
#' @param msg The message to echo
#'
#' @return a list
#' @export
#'
#' @examples
#' echo_back(msg = "I am a message and no one will stop me!")
echo_back <- function(msg = "") {
  list(msg = paste0("The message is: '", msg, "'"))
}


#' Plot a histogram of random numbers
#' @importFrom stats rnorm
#' @importFrom graphics hist
#' @return an image of a histogram
#' @export
#'
#' @examples
#' plot_histogram()
plot_histogram <- function() {
  rand <- rnorm(100)
  hist(rand)
}


#' Sum two numbers
#' @param a The first number to add
#' @param b The second number to add
#'
#' @return a number
#' @export
#'
#' @examples
#' sum_numbers(42, 114514)
sum_numbers <- function(a, b) {
  as.numeric(a) + as.numeric(b)
}
