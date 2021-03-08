library(iris3api)

#* @apiTitle Plumber Example API

#* Log requests
#* @filter log_request
iris3api::log_request

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
iris3api::echo_back

#* Plot a histogram
#* @png
#* @get /plot
iris3api::plot_histogram

#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
iris3api::sum_numbers
