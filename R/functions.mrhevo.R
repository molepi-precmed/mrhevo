#' Report a message to the terminal.
#'
#' @param mode One of \sQuote{info} (normal messages), \sQuote{note} (messages
#'        that require some highlighting), \sQuote{warn} (important information
#'        the user should definitely notice).
#' @param ... Strings to be reported.
#' @param LF Whether a newline character should be added at the end of the
#'        message (\code{TRUE} by default).
#'
#' @import crayon
#' @export
msg <- function(mode, ..., LF=TRUE) {
    message(mode(...), appendLF=LF)
}
info <- crayon::reset
note <- crayon::green
warn <- crayon::yellow
bold <- crayon::bold

#' Output a consisted header message.
#'
#' @param message Message to output.
#'
#' @export
header <- function(message) {
    msg(bold, note("\n#"), message)
}

#' Return upper bound on p-value where z is so extreme that pnorm(-z) returns 0.
#'
#' @param z Standard normal deviate.
#'
#' @return A string in scientific notation.
pnorm.extreme <- function(z, upper=TRUE) {
    ## https://www.johndcook.com/blog/norm-dist-bounds/
    ## get either upper bound on log10 tail prob or the lower bound
    c <- ifelse(upper, 8 / pi, 4)

    x <-  2 * sqrt(2 / pi) / (z + sqrt(z^2 + c))
    ln_p = -0.5 * z^2 + log(x)
    log10p <- ln_p / log(10)
    exponent <- floor(log10p)
    coeff <- 10^(log10p - exponent)
    string <- paste0(round(coeff), "E", exponent)
    return(string)
}

#' Reformat p-values that are in scientific notation for LaTeX.
#'
#' @param x Character vector of numbers in scientific notation.
#'
#' @return A character vector of LaTeX math expressions.
format.scinot.pvalue <- function(x) {
    x <- toupper(x)
    ## split x.split at E
    x.split <- as.numeric(unlist(strsplit(as.character(x), "E")))
    ## 1 significant figure
    x.split <- signif(as.numeric(x.split), 1)
    x.split <- t(matrix(x.split, nrow=2))
    ## handle cases where mantissa is rounded up to 10
    roundedto10 <- x.split[, 1]==10
    x.split[, 1][roundedto10] <- 1
    x.split[, 2][roundedto10] <- x.split[, 2][roundedto10] - 1
    p.latex <- sprintf("\\ensuremath{%.*f \\times 10^{%0*d}}",
                     0, x.split[, 1], 0, x.split[, 2])
    return(p.latex)
}

#' Generate formatted p-values from z values.
#'
#' @param z Standard normal deviate.
#' @param sigfig Number of significant figures to show in p-values.
#' @param neglogp.threshold.scinot Minus log10(pvalue) threshold for scientific
#'                                 notation.
#' @param neglogp.threshold Minus log10(pvalue) cutoff for thresholding
#'                          p-values.
#'
#' @return Vector of formatted p-values.
format.z.aspvalue <- function(z, sigfig=1, neglogp.threshold.scinot=3,
                              neglogp.threshold=NULL) {
    p <- signif(2 * pnorm(-abs(z)), sigfig)
    p.char <- toupper(as.character(p))
    ## pnorm.extreme returns a character string of form "NE-NNN"
    p.char[!is.na(p.char) & p.char=="0"] <-
        pnorm.extreme(z[!is.na(p.char) & p.char=="0"]) # where R outputs 0

    ## revert from scientific notation p-values above threshold for
    ## neglogp.threshold.scinot
    sci.revert <- grepl("E", p.char) & p > 10^-neglogp.threshold.scinot
    p.char[sci.revert] <-  format(p[sci.revert], scientific=FALSE)

    if (!is.null(neglogp.threshold)) { # thresholding of p values
        p.char[p < 10^-neglogp.threshold] <-
            paste0("<", format(10^-neglogp.threshold, scientific=FALSE))
    }
    ## format values in scientific notation for LaTeX
    p.char[grep("E", p.char)] <- format.scinot.pvalue(p.char[grep("E", p.char)])
    return(p.char)
}

#' Calculate summary stats for second step of two-step Mendelian randomization.
#'
#' @param Y Integer vector of binary outcome or numeric vector of continuous
#'          outcome.
#' @param Z Data.table of genetic instruments.
#' @param X_u Data.table of covariates.
#'
#' @return Data.table of coefficients for each genetic instrument
#'
#' @export
get_summarystatsforMR <- function(Y, Z, X_u) {
    # registerDoParallel(cores=10)
    YXZ.dt <- data.table(y=Y, Z, X_u)
    ## loop over instruments to fit regression of Y on Z, adjusted for X_u
    ## FIXME: implement a score test or parallelize
    coeffs <- foreach(i = 1:ncol(Z),
                      .combine=function(...) rbind(..., fill=TRUE),
                      .multicombine=TRUE) %dopar% {

        scoreid <- colnames(Z)[i]
        formula.string <- paste0("y ~ ", paste(colnames(X_u), collapse=" + "),
                                 " + ", scoreid)
        coeff <- summary(glm(data=YXZ.dt,
                             formula=as.formula(formula.string),
                             family="binomial"))$coefficients
        coeff <- as.data.table(coeff, keep.rownames="variable")
        coeff <- coeff[variable==scoreid]

    }
    colnames(coeffs) <- c("scoreid", "gamma_hat", "se.gamma_hat", "z", "p")
    return(coeffs)
}

#' Calculate ratios of coefficients.
#'
#' @param coeffs.dt Data table with coefficients for each genetic instrument.
#'
#' @return A data.table with coefficients for genetic instruments plus a column
#'         containing ratios of coefficients for each genetic instrument.
#'
#' @export
get_coeffratios <- function(coeffs.dt, use.delta=FALSE) {
    if (use.delta) {
        ## second-order Taylor expansions (delta method)
        coeffs.dt[, theta_IV :=
            gamma_hat / alpha_hat + se.alpha_hat^2 * gamma_hat / alpha_hat^3]
        coeffs.dt[, se.theta_IV :=
            sqrt((se.gamma_hat / alpha_hat)^2 +
                  gamma_hat^2 * se.alpha_hat^2 / alpha_hat^4)]
    } else {
        coeffs.dt[, theta_IV := gamma_hat / alpha_hat]  # ratio estimates
        coeffs.dt[, se.theta_IV := se.gamma_hat / alpha_hat]
    }
    coeffs.dt[, size.theta_IV := 0.3 * sum(se.theta_IV) / se.theta_IV]
    coeffs.dt[, inv.var := se.theta_IV^-2]
    return(coeffs.dt)
}

#' Calculate inverse-variance weighted MR estimator from summary stats.
#'
#' @param coeffs.dt Data.table with columns for coefficient and SE of
#'                  regression of exposure on instrument: \code{alpha_hat},
#'                  \code{se.alpha_hat}),
#'                  regression of outcome on instrument: \code{gamma_hat},
#'                  \code{se.gamma_hat} and coefficient ratios:
#'                  \code{theta_IV}, \code{se.theta_IV}.
#' @return A data.table with inverse-variance weighted estimate.
#'
#' @export
get_estimatorsMR <- function(coeffs.dt) {
    theta_IVW <- coeffs.dt[, sum(theta_IV * inv.var) / sum(inv.var)]
    se.theta_IVW  <- coeffs.dt[, sqrt(1 / sum(inv.var))]

    estimators.dt <- data.table(Estimator="Inverse variance weighted",
                                Estimate=theta_IVW,
                                SE=se.theta_IVW)
    estimators.dt[, z := Estimate / SE]
    estimators.dt[, pvalue := 2 * pnorm(-abs(z))]
    estimators.dt[, pvalue.formatted := format.z.aspvalue(z)]
    return(estimators.dt)
}

#' Calculate maximum likelihood estimate and p-value from posterior samples and
#' prior.
#'
#' @param x Numeric vector of samples from the posterior distribution of the
#'          parameter.
#' @param prior Numeric vector of values of the prior density for each element
#'              of x.
#' @param return.asplot Logical value determines whether the function returns
#'                      maximum likelihood estimate or a plot of the posterior
#'                      density and log-likelihood.
#'
#' @return If \code{return.asplot} is FALSE, return data.table with one row
#'         containing maximum likelihood estimate, standard error,
#'         test statistic, p-value, formatted p-value. Otherwise, create a plot
#'         of posterior density and log-likelihood.
#'
#' @import cowplot car
#' @export
mle.se.pval <- function(x, prior, return.asplot=FALSE) {
    invprior <- 1 / prior
    invprior <- invprior / sum(invprior)

    ## likelihood is posterior density divided by prior
    ## equivalently we can weight posterior samples by inverse of the prior
    ## when fitting a kernel density (usually Sheather-Jones is preferred to
    ## default bw)
    ## wider bandwidth gives better approximation to quadratic
    lik <- density(x, bw="SJ", adjust=2, weights=invprior)
    nonzero.lik <- which(lik$y > 0)
    logl <- log(lik$y[nonzero.lik])
    xvals <- lik$x[nonzero.lik]
    xvals.sq <- lik$x[nonzero.lik]^2

    ## possible refinement would be to weight the regression that is used to fit
    ## the quadratic approximation:
    fit.quad <- lm(logl ~ xvals + xvals.sq)
    a <- -as.numeric(fit.quad$coefficients[3]) # y = -a * x^2 + bx
    b <- as.numeric(fit.quad$coefficients[2]) # mle b/2a, se sqrt(1/2a)
    mle <- 0.5 * b / a
    stderr <- sqrt(0.5 / a)
    z <- mle / stderr
    pvalue <- 2 * pnorm(-abs(z))
    pvalue.formatted <- format.z.aspvalue(z)

    if (return.asplot) {
        loglik.dt <- data.table(logl.fit=logl, x=xvals,
                                logl.quad=-a * xvals^2 + b * xvals)
        loglik.dt[, rownum := .I]
        loglik.dt[, logl.fit := logl.fit - max(logl.fit)]
        loglik.dt[, logl.quad := logl.quad - max(logl.quad)]
        loglik.long <- melt(loglik.dt, id.vars="rownum",
                            measure.vars=c("logl.fit", "logl.quad"),
                            variable.name="curve", value.name="loglik")

        loglik.long <- loglik.dt[, .(rownum, x)][loglik.long, on="rownum"]
        loglik.long[, curve := car::recode(curve,
            "'logl.fit'='Smoothed curve'; 'logl.quad'='Fitted quadratic'",
            as.factor=TRUE)]

        p.loglik <- ggplot(loglik.long, aes(x=x, y=loglik, color=curve)) +
            geom_line() +
            xlab("Parameter value") +
            ylab("Log-likelihood") +
            theme(legend.position=c(0.5, 0.5)) +
            theme(legend.title=element_blank())

        xlimits <- c(min(loglik.long$x), max(loglik.long$x))
        fitted <- density(x)
        fitted.dt <- data.table(x=fitted$x, posterior=fitted$y)

        p.density <- ggplot(fitted.dt, aes(x=x, y=posterior)) +
            geom_line() +
            xlab("Parameter value") +
            ylab("Smoothed posterior density") +
            scale_x_continuous(limits=xlimits)

        p.bayesloglik <- cowplot::plot_grid(p.density, p.loglik, nrow=2)
        return(p.bayesloglik)
    } else {
        return(data.table(Estimate=mle, SE=stderr, z=z,
                          pvalue=pvalue, pvalue.formatted=pvalue.formatted))
    }
}

#' Run Stan model for Mendelian randomization with regularized horsehoe prior
#' on pleiotropic effects.
#'
#' @param sampling Logical. If set to \code{TRUE}, uses sampling rather than
#'                 variational approximation.
#' @param logistic Logical. If set to \code{TRUE}, uses logistic rather than
#'                 linear regression.
#' @param Z Data.table of genetic instruments.
#' @param Y Vector of outcomes.
#' @param sigma_y Standard deviation of outcome variable. Used only if
#'        \code{logistic=FALSE}.
#' @param X_u Data.table of unpenalized covariates.
#' @param alpha_hat Vector of estimated coefficients for effect of instruments
#'        on exposure.
#' @param se.alpha_hat Vector of standard errors for coefficients alpha_hat.
#' @param fraction_pleio Prior guess at fraction of instruments that have
#'        pleiotropic effects: values between 0.05 and 0.95 are allowed.
#' @param priorsd_theta Standard deviation of prior on theta.
#' @param model.dir Full path to STAN model directory.
#'
#' @returns An object of class stanfit.
#'
#' @example man/examples/runmrhevo.R
#' @export
#' @import data.table rstan ggplot2 bayesplot
run_mrhevo <- function(use.sampling=TRUE, logistic=TRUE,
                       Z, Y, sigma_y=1, X_u, alpha_hat, se.alpha_hat,
                       fraction_pleio=NULL, slab_scale=0.25, priorsd_theta=1,
                       vb.algo="meanfield", model.dir) {
    rstan_options(auto_write = TRUE)

    msg(bold, "Compiling stan model ... ")
    mr.stanmodel <- stan_model(file=file.path(model.dir, "MRHevo_logistic.stan"),
                               model_name="MRHevo.logistic", verbose=FALSE)
    msg(note, "Done.\n")

    ## check arguments for consistency
    stopifnot(length(unique(c(nrow(Z), nrow(X_u), length(Y)))) == 1)
    stopifnot(length(unique(c(ncol(Z), length(alpha_hat), length(se.alpha_hat)))) == 1)
    stopifnot(fraction_pleio >= 0.05 & fraction_pleio <= 0.95)

    N <- nrow(Z)
    J <- ncol(Z)
    X_u <- scale(X_u, center=TRUE, scale=TRUE)
    Z <- scale(Z, center=TRUE, scale=FALSE)

    ## priors
    ## weak prior on Y intercept
    scale_intercept_y <- 10

    ## prior sd of coeffs for unpenalized covariates X_u
    scale_beta_u <- 1

    ## prior scale of c is slab_scale

    ## 1 for half-Cauchy, large value specifies a gaussian prior on slab
    ## component
    slab_df <- 2

    ## 1 for half-Cauchy prior on global scale param: specifying a larger value
    ## will limit the narrowness of the spike component
    nu_global <- 1

    ## 1 for half-Cauchy, horseshoe+ or horseshoe if c is large
    nu_local <- 1

    ## prior guess of number of instruments that are pleiotropic
    if (is.null(fraction_pleio)) {
        fraction_pleio <- 0.5
    }
    r_pleio <- fraction_pleio * J

    ## Piironen and Vehtari recommend that the prior on tau should be chosen to
    ## have most of the prior mass
    ## near (r_pleio / (N * (J - r_pleio)) * sqrt(pseudovariance) / sqrt(N),
    ## where r_pleio is a prior guess for the number of nonzero coefficients

    ## prior median of a half-t distribution with nu_global df
    ## 0.82 with nu_global=1
    priormedian <- qt(p=0.75, df=nu_global, lower.tail=TRUE)

    if (logistic) {
        mu <- mean(Y)
        pseudovariance <- (mu * (1 - mu))^-1
        tau_0 <- (r_pleio / (J - r_pleio)) * sqrt(pseudovariance) / sqrt(N)
    } else {
        tau_0 <- (r_pleio / (J - r_pleio)) * sigma_y / sqrt(N)
    }

    ## choose prior on tau so that most of the prior mass is near tau_0
    ## choose scale_global so that the prior median of the half-t distribution
    ## equates to tau_0
    scale_global <- tau_0 / priormedian

    ## sample the posterior
    data.stan <- list(logistic=as.integer(logistic),
                    Z=as.matrix(Z), Y=Y, X_u=as.matrix(X_u),
                    N=N, J=J, U=ncol(X_u),
                    alpha_hat=alpha_hat,
                    sd_alpha_hat=se.alpha_hat,
                    nu_global=nu_global,
                    nu_local=nu_local,
                    scale_global=scale_global,
                    scale_intercept_y=scale_intercept_y,
                    scale_beta_u=scale_beta_u,
                    priorsd_theta=priorsd_theta,
                    slab_scale=slab_scale,
                    slab_df=slab_df, priorsd_theta=priorsd_theta)

    if (use.sampling) {
        fit.mc <- rstan::sampling(object=mr.stanmodel,
                                  data=data.stan,
                                  iter=1200, warmup=400,
                                  cores=4,
                                  chains=4,
                                  refresh=200,
                                  control=list(adapt_delta=0.95),
                                  verbose=FALSE)
    } else {
        fit.mc <- rstan::vb(object=mr.stanmodel,
                            data=data.stan,
                            algorithm=vb.algo,
                            refresh=5000,
                            iter=20000,
                            adapt_engaged=TRUE,
                            tol_rel_obj=0.01)
    }
    return(fit.mc)
}

#' Run Stan model for Mendelian randomization with regularized horsehoe prior
#' on pleiotropic effects, using summary statistics only.
#'
#' @param alpha_hat Vector of estimated coefficients for effect of instruments
#'        on exposure.
#' @param se.alpha_hat Vector of standard errors for coefficients alpha_hat.
#' @param gamma_hat Vector of estimated coefficients for effect of instruments
#'        on outcome.
#' @param se.gamma_hat Vector of standard errors for coefficients gamma_hat.
#' @param fraction_pleio Prior guess at fraction of instruments that have
#'        pleiotropic effects: values between 0.05 and 0.95 are allowed.
#' @param slab_scale scale param of prior on direct effects.
#' @param priorsd_theta Standard deviation of prior on theta.
#' @param model.dir Full path to STAN model directory.
#'
#' @return An object of class stanfit.
#' @export
run_mrhevo.sstats <- function(alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat,
                              fraction_pleio=NULL, slab_scale=0.2, slab_df=2,
                              priorsd_theta=1, model.dir,
                              hierarchical_alpha = TRUE) {
    require(rstan)
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)

    msg(bold, "Compiling stan model ... ")
    mr.sstats.stanmodel <- stan_model(file=file.path(model.dir,
                                                     "MRHevo_summarystats.stan"),
                                      model_name="MRHevo.summarystats",
                                      verbose=FALSE)
    msg(note, "Done.\n")

    ## check arguments for consistency
    stopifnot(length(alpha_hat)==length(gamma_hat))
    stopifnot(fraction_pleio >= 0.01 & fraction_pleio <= 0.99)

    J <- length(alpha_hat)

    ## Taylor expansion
    var.theta_IV_delta <- (se.gamma_hat / alpha_hat)^2 +
    gamma_hat^2 * se.alpha_hat^2 / alpha_hat^4

    ## mean Fisher info for IV estimate
    info <- mean( 1 / var.theta_IV_delta)

    ## df of half-t priors
    ## 1 for half-Cauchy prior on global scale param: specifying a larger value
    ## will limit the narrowness of the spike component
    nu_global <- 1
    ## 1 for half-Cauchy, horseshoe+ or horseshoe if c is large
    nu_local <- 1

    ## prior median of a half-t distribution with nu_global df
    tau0 <- set.tau0(fraction_pleio=fraction_pleio, nu_global=nu_global, J=J,
                     info=info)

    ## choose prior on tau so that most of the prior mass is near tau_0
    ## choose scale_global so that the prior median of the half-t distribution
    ## equates to tau_0
    ## 0.82 with nu_global=1
    priormedian <- qt(p=0.75, df=nu_global, lower.tail=TRUE)
    scale_global <- tau0 / priormedian

    data.stan <- list(J=length(alpha_hat),
                      gamma_hat=gamma_hat,
                      sd_gamma_hat=se.gamma_hat,
                      alpha_hat=alpha_hat,
                      sd_alpha_hat=se.alpha_hat,
                      nu_global=nu_global,
                      nu_local=nu_local,
                      scale_global=scale_global,
                      priorsd_theta=priorsd_theta,
                      slab_scale=slab_scale,
                      slab_df=slab_df,
                      hierarchical_alpha=as.integer(hierarchical_alpha))
    msg(bold, "Sampling posterior distribution ... ")
    fit.mc <- rstan::sampling(object=mr.sstats.stanmodel,
                              data=data.stan,
                              iter=3000, warmup=1000,
                              cores=4,
                              chains=4,
                              refresh=1000,
                              control=list(adapt_delta=0.99),
                              verbose=TRUE)
    msg(note, "Done.\n")

    ## Check for sampling warnings
    sampler_params <- rstan::get_sampler_params(fit.mc, inc_warmup = FALSE)
    divergent <- sum(sapply(sampler_params, function(x) sum(x[,"divergent__"])))
    total_iterations <- sum(sapply(sampler_params, nrow))

    if (divergent > 0) {
      warning(paste0("There were ", divergent, " divergent transitions after warmup (",
                     round(100 * divergent / total_iterations, 1), "% of iterations). ",
                     "This may indicate sampling difficulties. To address this:\n",
                     "  - Reduce slab_scale (try 0.1 or 0.05)\n",
                     "  - Increase slab_df (try 4 or 8)\n",
                     "  - Increase warmup iterations (e.g., warmup=2000)"))
    }

    return(fit.mc)
}

#' Run Stan model for Mendelian randomization with regularized horsehoe prior
#' on pleiotropic effects, using summary statistics only.
#'
#' @param alpha_hat Vector of estimated coefficients for effect of instruments
#'        on exposure.
#' @param se.alpha_hat Vector of standard errors for coefficients alpha_hat.
#' @param gamma_hat Vector of estimated coefficients for effect of instruments
#'        on outcome.
#' @param se.gamma_hat Vector of standard errors for coefficients gamma_hat.
#' @param fraction_pleio Prior guess at fraction of instruments that have
#'        pleiotropic effects: values between 0.05 and 0.95 are allowed.
#' @param slab_scale scale param of prior on direct effects.
#' @param priorsd_theta Standard deviation of prior on theta.
#' @param model.dir Full path to STAN model directory.
#'
#' @return An object of class stanfit.
#' @export
run_mrhevo.fixedtau <- function(alpha_hat, se.alpha_hat, gamma_hat,
                                se.gamma_hat, tau=1E-6, slab_scale=0.2,
                                slab_df=2, priorsd_theta=1, model.dir) {
    rstan_options(auto_write = TRUE)

    ## check arguments for consistency
    stopifnot(length(alpha_hat)==length(gamma_hat))

    msg(bold, "Compiling stan model ... ")
    mr.fixedtau.stanmodel <- stan_model(file=file.path(model.dir,
                                                       "MRHevo_fixedtau.stan"),
                                        model_name="MRHevo.fixedtau",
                                        verbose=FALSE)
    msg(note, "Done.\n")

    J <- length(alpha_hat)

    ## priors
    ## 1 for half-Cauchy, horseshoe+ or horseshoe if c is large
    nu_local <- 1

    data.stan <- list(J=length(alpha_hat),
                      gamma_hat=gamma_hat,
                      sd_gamma_hat=se.gamma_hat,
                      alpha_hat=alpha_hat,
                      sd_alpha_hat=se.alpha_hat,
                      tau=tau,
                      nu_local=nu_local,
                      priorsd_theta=priorsd_theta,
                      slab_scale=slab_scale,
                      slab_df=slab_df, priorsd_theta=priorsd_theta)

    msg(bold, "Sampling posterior distribution ... ")
    fit.mc <- rstan::sampling(object=mr.fixedtau.stanmodel,
                              data=data.stan,
                              iter=3000, warmup=1000,
                              cores=4,
                              chains=4,
                              refresh=1000,
                              control=list(adapt_delta=0.99),
                              verbose=TRUE)
    msg(note, "Done.\n")
    return(fit.mc)
}

#' Get value \code{tau0} for global shrinkage parameter \code{tau} given
#' expectation of \code{fraction_pleio}.
#'
#' @param fraction_pleio Prior guess at fraction of effects that are nonzero.
#' @param nu_global Shape parameter of gamma prior.
#' @param J number of variables.
#' @param pseudovariance Variance of outcome variable.
#'
#' @return Value of tau0.
set.tau0 <- function(fraction_pleio=NULL, nu_global=1, J, info) {
    ## prior guess of number of instruments that are pleiotropic
    if (is.null(fraction_pleio)) {
        fraction_pleio <- 0.5
    }
    r_pleio <- fraction_pleio * J

    ## Piironen and Vehtari recommend that the prior on tau should be chosen to
    ## have most of the prior mass
    ## near (r_pleio / (N * (J - r_pleio)) * sqrt(pseudovariance) / sqrt(N),
    ## where r_pleio is a prior guess for the number of nonzero coefficients
    tau0 <- (r_pleio / (J - r_pleio)) / sqrt(info)
    return(tau0)
}

#' Plot instrumental variable estimates with MLE as slope line.
#'
#' @param alpha_hat Vector of estimated coefficients for effect of instruments
#'        on exposure.
#' @param se.alpha_hat Vector of standard errors for coefficients alpha_hat.
#' @param gamma_hat Vector of estimated coefficients for effect of instruments
#'        on outcome.
#' @param se.gamma_hat Vector of standard errors for coefficients gamma_hat.
#' @param theta Maximum likelihood estimate of causal effect (slope).
#' @param qtlname Optional vector of QTL names for labeling points.
#'
#' @return A ggplot object showing IV estimates with MLE as line through origin.
#'
#' @import ggplot2 data.table ggrepel
#' @export
plot_iv_estimates <- function(alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat, theta, qtlname = NULL) {
    theta_IV <- gamma_hat / alpha_hat
    se.theta_IV <- se.gamma_hat / abs(alpha_hat)

    iv.dt <- data.table(
        alpha_hat = alpha_hat,
        gamma_hat = gamma_hat,
        theta_IV = theta_IV,
        se.theta_IV = se.theta_IV
    )

    if (!is.null(qtlname)) {
        iv.dt[, qtlname := qtlname]
    }

    p <- ggplot(iv.dt, aes(x = alpha_hat, y = gamma_hat)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = theta,
                   color = "red", linewidth = 0.8, linetype = "dashed") +
        xlab("Effect of instrument on exposure (alpha)") +
        ylab("Effect of instrument on outcome (gamma)") +
        scale_x_continuous(limits = c(0, NA)) +
        theme_bw() +
        theme(legend.position = "none")

    if (!is.null(qtlname)) {
        p <- p + geom_text_repel(aes(label = qtlname), size = 2.5, max.overlaps = 20)
    }

    return(p)
}

#' Plot pairs of posterior samples.
#'
#' @param fit Stan fit object or mrhevo_numpyro object.
#' @param pars Vector of parameter names to plot (default c("theta", "f", "log_tau")).
#'
#' @return A ggplot object showing pairs plot.
#'
#' @import bayesplot ggplot2 cowplot
#' @export
plot_posterior_pairs <- function(fit, pars = c("theta", "f", "log_tau")) {
    if (inherits(fit, "mrhevo_numpyro")) {
        # NumPyro output - extract directly from posterior list
        post <- fit$posterior
        
        # Build data frame for requested parameters
        df_list <- list()
        
        for (p in pars) {
            if (p %in% names(post)) {
                vals <- post[[p]]
                # Flatten to vector
                df_list[[p]] <- as.vector(vals)
            }
        }
        
        df_plot <- as.data.frame(df_list)
        
        if (ncol(df_plot) < 2) {
            stop("Need at least 2 parameters for pairs plot")
        }
        
        # Create pairs plot using cowplot grid
        plots_list <- list()
        n_pars <- ncol(df_plot)
        
        for (i in 1:n_pars) {
            for (j in 1:n_pars) {
                if (i == j) {
                    # Diagonal: histogram
                    p_ij <- ggplot(df_plot, aes(x = .data[[names(df_plot)[i]]])) +
                        geom_histogram(fill = "steelblue", color = "white", bins = 30) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(text = ggplot2::element_text(size = 8))
                } else if (j < i) {
                    # Lower triangle: scatter plot
                    p_ij <- ggplot(df_plot, aes(x = .data[[names(df_plot)[j]]], y = .data[[names(df_plot)[i]]])) +
                        geom_point(alpha = 0.2, size = 0.3) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(text = ggplot2::element_text(size = 8))
                } else {
                    # Upper triangle: correlation
                    cor_val <- cor(df_plot[, j], df_plot[, i], use = "complete.obs")
                    p_ij <- ggplot() + 
                        annotate("text", x = 0.5, y = 0.5, label = sprintf("r = %.2f", cor_val), size = 3) +
                        theme_void() +
                        ggplot2::theme(text = ggplot2::element_text(size = 8))
                }
                plots_list[[paste0(i, "_", j)]] <- p_ij
            }
        }
        
        # Arrange in grid
        p <- cowplot::plot_grid(plotlist = plots_list, 
                                 nrow = n_pars, ncol = n_pars,
                                 labels = names(df_plot),
                                 label_size = 8)
        
        return(p)
    } else {
        # Stan output - use bayesplot
        np <- bayesplot::nuts_params(fit, type = "divergent")
        
        p <- bayesplot::mcmc_pairs(fit, pars = pars,
                                   np = np,
                                   off_diag_args = list(alpha = 0.3))
        
        p <- p + ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 18))
        
        return(p)
    }
}

#' Plot histogram of kappa shrinkage coefficients.
#'
#' @param fit Stan fit object or mrhevo_numpyro object.
#' @param bin_width Width of bins for histogram (default 0.02).
#'
#' @return A ggplot object showing histogram of kappa values.
#'
#' @import ggplot2
#' @export
plot_kappa_hist <- function(fit, bin_width = 0.02) {
    if (inherits(fit, "mrhevo_numpyro")) {
        # NumPyro output
        kappa <- fit$posterior$kappa
    } else {
        # Stan output
        posterior <- rstan::extract(fit)
        kappa <- posterior$kappa
    }

    kappa_vec <- as.vector(kappa)
    kappa_dt <- data.table::data.table(kappa = kappa_vec)

    p <- ggplot(kappa_dt, aes(x = kappa, y = after_stat(density))) +
        geom_histogram(binwidth = bin_width, fill = "steelblue", color = "white",
                       boundary = 0) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        xlab("Kappa (shrinkage coefficient)") +
        ylab("Proportion") +
        ggtitle(paste("Histogram of kappa (bin width =", bin_width, ")")) +
        theme_bw() +
        theme(legend.position = "none")

    return(p)
}

#' Run Stan model for Mendelian randomization using NumPyro.
#'
#' @param alpha_hat Vector of estimated coefficients for effect of instruments
#'        on exposure.
#' @param se.alpha_hat Vector of standard errors for coefficients alpha_hat.
#' @param gamma_hat Vector of estimated coefficients for effect of instruments
#'        on outcome.
#' @param se.gamma_hat Vector of standard errors for coefficients gamma_hat.
#' @param fraction_pleio Prior guess at fraction of instruments that have
#'        pleiotropic effects.
#' @param slab_scale Scale for the slab component of regularized horseshoe.
#' @param slab_df Degrees of freedom for slab component.
#' @param priorsd_theta Standard deviation of prior on theta.
#' @param model_path Path to the mrhevo_pyro.py script.
#' @param env_path Path to Python virtual environment (default: ~/.virtualenvs/mrhevo).
#' @param num_warmup Number of warmup iterations (default 500).
#' @param num_samples Number of samples (default 1000).
#' @param num_chains Number of chains (default 4).
#' @param target_accept_prob Target acceptance probability (default 0.95).
#'
#' @return A list with posterior samples in rstan-like format.
#'
#' @import reticulate
#' @export
run_mrhevo.numpyro <- function(alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat,
                                fraction_pleio = 0.5, slab_scale = 0.2, slab_df = 2,
                                priorsd_theta = 1, model_path = NULL,
                                env_path = NULL,
                                hierarchical_alpha = TRUE,
                                num_warmup = 500, num_samples = 1000,
                                num_chains = 4, target_accept_prob = 0.95) {

    if (is.null(env_path)) {
        env_path <- "~/.virtualenvs/mrhevo"
    }

    reticulate::use_virtualenv(env_path, required = TRUE)

    if (is.null(model_path)) {
        model_path <- system.file("python", "mrhevo_pyro.py", package = "mrhevo")
        if (model_path == "") {
            stop("mrhevo_pyro.py not found. Please provide model_path or install the Python script.")
        }
    }

    msg(bold, "Setting up NumPyro...")

    np <- reticulate::import("numpy", as = "np")
    jnp <- reticulate::import("jax.numpy", as = "jnp")
    random <- reticulate::import("jax.random")
    numpyro <- reticulate::import("numpyro")
    numpyro_distributions <- reticulate::import("numpyro.distributions")

    reticulate::source_python(model_path)

    J <- length(alpha_hat)
    var.theta_IV_delta <- (se.gamma_hat / alpha_hat)^2 + gamma_hat^2 * se.alpha_hat^2 / alpha_hat^4
    info <- mean(1 / var.theta_IV_delta)

    nu_global <- 1
    nu_local <- 1
    r_pleio <- fraction_pleio * J
    tau0 <- (r_pleio / (J - r_pleio)) / sqrt(info)
    priormedian <- qt(p = 0.75, df = nu_global, lower.tail = TRUE)
    scale_global <- tau0 / priormedian

    msg(bold, "Running NumPyro MCMC...")

    model <- mrhevo(
        alpha_hat = np$array(alpha_hat),
        se_alpha_hat = np$array(se.alpha_hat),
        gamma_hat = np$array(gamma_hat),
        se_gamma_hat = np$array(se.gamma_hat),
        slab_scale = slab_scale,
        slab_df = slab_df,
        scale_global = scale_global,
        priorsd_theta = priorsd_theta,
        hierarchical_alpha = as.integer(hierarchical_alpha)
    )

    nuts_kernel <- numpyro$infer$MCMC(
        numpyro$infer$NUTS(model, target_accept_prob = target_accept_prob),
        num_warmup = as.integer(num_warmup),
        num_samples = as.integer(num_samples),
        num_chains = as.integer(num_chains),
        chain_method = "parallel",
        progress_bar = TRUE
    )

    start_time <- Sys.time()
    rng_key <- random$PRNGKey(42L)
    nuts_kernel$run(rng_key,
                    alpha_hat = np$array(alpha_hat),
                    se_alpha_hat = np$array(se.alpha_hat),
                    gamma_hat = np$array(gamma_hat),
                    se_gamma_hat = np$array(se.gamma_hat),
                    slab_scale = as.double(slab_scale),
                    slab_df = as.double(slab_df),
                    scale_global = as.double(scale_global),
                    priorsd_theta = as.double(priorsd_theta),
                    hierarchical_alpha = as.integer(hierarchical_alpha))
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    msg(note, paste0("NumPyro sampling completed in ", round(elapsed, 1), " seconds\n"))

    samples <- nuts_kernel$get_samples()

    sample_names <- c("theta", "tau", "log_tau", "eta", "log_eta", "f", "b",
                      "alpha", "beta", "kappa", "lambda_tilde")

    np <- reticulate::import("numpy", as = "np")
    az <- reticulate::import("arviz")

    idata <- az$from_numpyro(nuts_kernel)

    posterior_list <- list()
    for (name in sample_names) {
        tryCatch({
            if (name %in% names(idata$posterior)) {
                arr <- np$asarray(idata$posterior[[name]])
                arr <- reticulate::py_to_r(arr)
                if (!is.null(dim(arr))) {
                    if (length(dim(arr)) == 2) {
                        arr <- as.vector(arr)
                    }
                } else {
                    arr <- as.vector(arr)
                }
                posterior_list[[name]] <- arr
            }
        }, error = function(e) {
            cat("Error extracting", name, ":", conditionMessage(e), "\n")
        })
    }

    fit_numpyro <- list(
        posterior = posterior_list,
        elapsed = elapsed,
        num_chains = num_chains,
        num_samples = num_samples,
        num_warmup = num_warmup
    )

    class(fit_numpyro) <- "mrhevo_numpyro"

    return(fit_numpyro)
}

#' Convert NumPyro posterior to rstan-like format.
#'
#' @param fit A mrhevo_numpyro object.
#'
#' @return A list with posterior samples in rstan-compatible format.
#'
#' @export
convert_to_rstan <- function(fit) {
    if (!inherits(fit, "mrhevo_numpyro")) {
        stop("fit must be a mrhevo_numpyro object")
    }

    J <- ncol(fit$posterior$kappa)
    n_samples <- dim(fit$posterior$theta)[1]
    n_chains <- dim(fit$posterior$theta)[2]

    dim(fit$posterior$theta) <- c(n_samples * n_chains)
    dim(fit$posterior$f) <- c(1)
    dim(fit$posterior$log_tau) <- c(1)
    dim(fit$posterior$log_eta) <- c(1)

    result <- list(
        theta = fit$posterior$theta,
        f = fit$posterior$f,
        log_c = fit$posterior$log_eta,
        log_tau = fit$posterior$log_tau,
        kappa = fit$posterior$kappa,
        alpha = fit$posterior$alpha,
        beta = fit$posterior$beta
    )

    return(result)
}
