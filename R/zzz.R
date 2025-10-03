# Package initialization

.onLoad <- function(libname, pkgname) {
  # This function ensures the package loads properly
  # It helps resolve issues where functions aren't immediately available

  # Force namespace evaluation to ensure all exports are available
  # This can help with lazy loading issues
  ns <- getNamespace(pkgname)

  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Optional startup message
  # packageStartupMessage("sir package (version ", utils::packageVersion(pkgname), ") loaded")
  invisible()
}