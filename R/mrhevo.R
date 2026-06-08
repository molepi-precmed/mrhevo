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

## Suppress R CMD check notes for data.table column names and stats functions
## used inside data.table expressions or ggplot2 aesthetics.
utils::globalVariables(c(
    ".", "i",
    "theta_IV", "gamma_hat", "alpha_hat", "se.alpha_hat", "se.gamma_hat",
    "se.theta_IV", "size.theta_IV", "inv.var",
    "Estimate", "SE", "z", "pvalue", "pvalue.formatted",
    "variable",
    "rownum", "logl.fit", "logl.quad", "loglik", "curve", "posterior",
    "pt_color"
))

#' @import doParallel
#' @import foreach
#' @importFrom utils packageVersion
#' @importFrom stats pnorm qt lm density glm as.formula cor
.onAttach <- function(libname, pkgname) {
    ## Do not spawn parallel workers during R CMD check (CRAN limits simultaneous processes).
    check_env <- identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")
    if (is.null(getOption("cores")))
        options(cores = if (check_env) 1L
                        else min(floor(parallel::detectCores() / 2), 10))
    if (!check_env)
        doParallel::registerDoParallel()

    packageStartupMessage("MRHEVO ", packageVersion("mrhevo"), ":")
    packageStartupMessage("    Currently using ", foreach::getDoParWorkers(),
                          " cores, set 'options(cores=<n.cores>)' to change.")
    packageStartupMessage("    NumPyro backend: call install_mrhevo_python() to set up.")
    packageStartupMessage("    Stan backend: call install_mrhevo_stan() to set up.")
}


#' Install Python dependencies for the NumPyro backend.
#'
#' Creates a Python virtual environment and installs \code{jax},
#' \code{numpyro}, and \code{arviz}.  Before installing, checks whether the
#' CPU supports AVX2 instructions required by JAX.  If the CPU is
#' incompatible, warns the user and offers to install the Stan backend instead
#' via \code{\link{install_mrhevo_stan}}.
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

    ## Check CPU compatibility: JAX requires AVX2 (x86_64 CPUs since ~2013).
    jax_ok <- .check_jax_cpu_compat()
    if (!jax_ok) {
        warning("This CPU does not support AVX2 instructions required by JAX/NumPyro.\n",
                "The NumPyro backend will not work on this machine.")
        if (interactive()) {
            answer <- readline("Install the Stan backend instead? [y/N] ")
            if (tolower(trimws(answer)) %in% c("y", "yes")) {
                install_mrhevo_stan()
                return(invisible(envpath))
            }
        }
        message("Run install_mrhevo_stan() to use the Stan backend instead.")
        return(invisible(envpath))
    }

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

## Internal helper: returns TRUE if CPU supports AVX2 (required by JAX).
.check_jax_cpu_compat <- function() {
    cpu_flags <- tryCatch({
        if (file.exists("/proc/cpuinfo")) {
            paste(readLines("/proc/cpuinfo", warn=FALSE), collapse=" ")
        } else {
            ## macOS fallback
            paste(system("sysctl -n machdep.cpu.features machdep.cpu.leaf7_features 2>/dev/null",
                         intern=TRUE, ignore.stderr=TRUE), collapse=" ")
        }
    }, error=function(e) "")
    grepl("avx2", tolower(cpu_flags))
}

#' Install the Stan backend for mrhevo.
#'
#' Installs the \pkg{rstan} and \pkg{bayesplot} R packages, which are only
#' required when using the Stan sampler (\code{\link{run_mrhevo_stan}}).
#' Use this function if the NumPyro backend is unavailable (e.g., because
#' the CPU does not support AVX2); see \code{\link{install_mrhevo_python}}.
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

#' Example dataset for MRHevo analysis.
#'
#' A data.table containing summary statistics for genetic instruments used
#' in an exemplar Mendelian randomisation analysis.
#'
#' @format A data.table with columns:
#' \describe{
#'   \item{scoreid}{Identifier for the genetic score/instrument.}
#'   \item{qtlname}{Name of the quantitative trait locus.}
#'   \item{alpha_hat}{Estimated coefficient for effect of instrument on exposure.}
#'   \item{se.alpha_hat}{Standard error of \code{alpha_hat}.}
#'   \item{gamma_hat}{Estimated coefficient for effect of instrument on outcome.}
#'   \item{se.gamma_hat}{Standard error of \code{gamma_hat}.}
#' }
#' @source Exemplar dataset bundled with the mrhevo package.
"coeffs.dt"
