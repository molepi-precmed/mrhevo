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
    packageStartupMessage("    NumPyro (default backend): call install_mrhevo_python() to set up.")
    packageStartupMessage("    Stan (optional backend): call install_mrhevo_stan() to set up.")
}

.onLoad <- function(libname, pkgname) {
    ## Install NumPyro Python environment on first load if not already present
    env_path <- path.expand("~/.virtualenvs/mrhevo")
    if (!reticulate::virtualenv_exists(env_path)) {
        tryCatch(
            install_mrhevo_python(envpath=env_path, ask=FALSE),
            error=function(e) {
                warning("Could not auto-install NumPyro environment: ",
                        conditionMessage(e),
                        "\nRun install_mrhevo_python() manually to set up.")
            }
        )
    }
}

#' Install Python dependencies for the NumPyro backend.
#'
#' Creates a virtual environment and installs \code{jax}, \code{numpyro},
#' and \code{arviz}.  This is called automatically on first package load;
#' call it manually if the environment needs to be rebuilt.
#'
#' @param envpath Path for the virtual environment
#'        (default: \code{~/.virtualenvs/mrhevo}).
#' @param ask Logical. If \code{TRUE} (the default when called interactively),
#'        prompt before creating the environment.
#'
#' @return Invisibly returns the environment path.
#'
#' @import reticulate
#' @export
install_mrhevo_python <- function(envpath="~/.virtualenvs/mrhevo", ask=TRUE) {
    envpath <- path.expand(envpath)
    if (ask && interactive()) {
        answer <- readline(
            paste0("Install NumPyro Python environment at ", envpath, "? [y/N] "))
        if (!tolower(trimws(answer)) %in% c("y", "yes")) {
            message("Skipping NumPyro installation.")
            return(invisible(envpath))
        }
    }
    message("Creating virtual environment at ", envpath, " ...")
    reticulate::virtualenv_create(envname=envpath)
    message("Installing jax, numpyro, arviz ...")
    reticulate::virtualenv_install(
        envname=envpath,
        packages=c("jax", "numpyro", "arviz"),
        ignore_installed=FALSE)
    message("NumPyro environment ready.")
    invisible(envpath)
}

#' Install the Stan backend for mrhevo.
#'
#' Installs the \pkg{rstan} and \pkg{bayesplot} R packages, which are only
#' required when using the Stan sampler (\code{\link{run_mrhevo.sstats}} with
#' \code{model.dir} pointing to the bundled Stan models).
#'
#' @param repos CRAN mirror to use (defaults to the current option).
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
install_mrhevo_stan <- function(repos=getOption("repos")) {
    pkgs <- c("rstan", "bayesplot")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
    if (length(missing) == 0) {
        message("Stan packages already installed: ", paste(pkgs, collapse=", "))
        return(invisible(NULL))
    }
    message("Installing: ", paste(missing, collapse=", "), " ...")
    utils::install.packages(missing, repos=repos)
    invisible(NULL)
}
