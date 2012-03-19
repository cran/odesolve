.onAttach <- function(libname, pkgname) {
  packageStartupMessage("odesolve is deprecated!  Use the solvers in deSolve instead.\nodesolve will be removed from CRAN by the end of 2012.")
}
