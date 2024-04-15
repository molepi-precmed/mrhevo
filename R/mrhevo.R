## mrhevo.R
##
## Package documentation and namespace hooks.
##

#' Inference of causal effects by Mendelian randomisation, marginalising over
#' distribution of direct effects
#'
#' @author
#' Paul McKeigue \email{paul.mckeigue@@ed.ac.uk},
#' Andrii Iakovliev \email{andrii.iakovliev@@ed.ac.uk},
#' Buddhi Erabadda \email{b.erabadda@@ed.ac.uk},
#' Athina Spiliopoulou \email{a.spiliopoulou@@ed.ac.uk}
#'
#' @docType package
"_PACKAGE"

## Environment for package global variables
mrhevo.env <- new.env(parent=emptyenv())

#' @import doParallel
#' @import foreach
#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    ## put a cap on the number of cores used by default
    if (is.null(getOption("cores")))
        options(cores=min(floor(parallel::detectCores() / 2), 10))
    doParallel::registerDoParallel()

    packageStartupMessage("MRHEVO ", packageVersion("mrhevo"), ":")
    packageStartupMessage("    Currently using ", foreach::getDoParWorkers(),
                          " cores, set 'options(cores=<n.cores>)' to change.")
}
