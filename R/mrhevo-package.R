#' The 'mrhevo' package.
#'
#' @description Inference of causal effects by Mendelian randomisation,
#'              marginalising over distribution of direct effects.
#'
#' @author
#' Paul McKeigue \email{paul.mckeigue@@ed.ac.uk},
#' Andrii Iakovliev \email{andrii.iakovliev@@ed.ac.uk},
#' Buddhi Erabadda \email{b.erabadda@@ed.ac.uk},
#' Athina Spiliopoulou \email{a.spiliopoulou@@ed.ac.uk}
#'
#' @name mrhevo-package
#' @aliases mrhevo
#' @useDynLib mrhevo, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import doParallel
#' @import foreach
#' @importFrom utils packageVersion
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.3. https://mc-stan.org
#'
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) Rcpp::loadModule(m, what = TRUE)
}

.onAttach <- function(libname, pkgname) {
    ## put a cap on the number of cores used by default
    if (is.null(getOption("cores")))
        options(cores=min(floor(parallel::detectCores() / 2), 10))
    doParallel::registerDoParallel()

    packageStartupMessage("MRHEVO ", packageVersion("mrhevo"), ":")
    packageStartupMessage("    Currently using ", foreach::getDoParWorkers(),
                          " cores, set 'options(cores=<n.cores>)' to change.")
}
