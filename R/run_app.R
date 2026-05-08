#' Launch the nestprepper Shiny App
#' 
#' @description This launches a web-based dashboard for QAQC and modeling.
#' @importFrom shiny runApp
#' @export
launch_app <- function() {
  app_dir <- system.file("shiny-examples", "nestprepper-app", package = "nestprepper")
  
  if (app_dir == "") {
    stop("Could not find the app directory. Try re-installing `nestprepper`.", call. = FALSE)
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
}