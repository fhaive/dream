#' Start the Flask API server
#'
#' This function starts a server using the specified Python interpreter.
#' It creates a new process and runs the "dream.api" module using the "-m" flag.
#' The standard output and standard error of the process are redirected to two log files.
#'
#' @param python Path to the Python interpreter (default is "python")
#' @return An instance of the "process" class representing the server process
#' @export
start_server <- function(python="python") {
  
  p <- processx::process$new(command=python, args=c("-m", "dream.api"), stdout="dream_api_out.log", stderr="dream_api_err.log")
  Sys.sleep(15)
  
  if(p$is_alive()) {
    message("server has been started")
  }
  
  return(p)
}

#' Stop the Flask API server
#'
#' This function stops the server process. It sends an interrupt signal to the
#' process to gracefully terminate it.
#'
#' @param p An instance of the "process" class representing the server process
#' @export
stop_server <- function(p) {
  if(p$is_alive()) {
    p$interrupt()
    message("server has been shutdown")
  }
  else {
    message("server process is not alive")
  }
}