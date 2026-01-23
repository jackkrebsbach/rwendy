.onLoad <- function(libname, pkgname) {

  # Set up package-specific logger namespace with INFO as default

  logger::log_threshold(logger::INFO, namespace = pkgname)
}

#' Set wendy logging level
#'
#' Control the verbosity of logging output from the wendy package.
#'
#' @param level Character string. One of "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "OFF".
#'   Default is "INFO".
#'
#' @return Invisible NULL. Called for side effect of setting log level.
#'
#' @examples
#' \dontrun{
#' # Enable debug output
#' wendy_log_level("DEBUG")
#'
#' # Silence all logging
#' wendy_log_level("OFF")
#'
#' # Reset to default
#' wendy_log_level("INFO")
#' }
#'
#' @export
wendy_log_level <- function(level = c("INFO", "DEBUG", "TRACE", "WARN", "ERROR", "OFF")) {
  level <- match.arg(level)
  logger::log_threshold(level, namespace = "wendy")
  invisible(NULL)
}
