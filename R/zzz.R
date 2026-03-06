# package initialization

.onLoad <- function(libname, pkgname) {
  # force namespace evaluation so all exports are available
  ns <- getNamespace(pkgname)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  invisible()
}