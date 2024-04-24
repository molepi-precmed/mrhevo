source("runmrhevo.R")


weighted.median.sd.bowden <- function(alpha_hat, gamma_hat, se.alpha_hat,
                                      se.gamma_hat, weights) {
    med = numeric(1000)
    for (i in 1:1000) {
        alpha_hat.boot = rnorm(length(alpha_hat), mean=alpha_hat,
                               sd=se.alpha_hat)
        gamma_hat.boot = rnorm(length(gamma_hat), mean=gamma_hat,
                               sd=se.gamma_hat)
        betaIV.boot = gamma_hat.boot / alpha_hat.boot
        med[i] = matrixStats::weightedMedian(betaIV.boot, weights)
    }
    return(sd(med))
}

use.delta <- TRUE

tau <-  mrhevo_samples.list[["tau"]]
lambda_tilde <- mrhevo_samples.list[["lambda_tilde"]] # numdraws x J
alpha.sim <- mrhevo_samples.list[["alpha"]] # numdraws x J

numdraws <- length(tau)
J <- ncol(lambda_tilde)

## simulate beta values from distribution
beta.sim <- matrix(rnorm(n=numdraws * J), nrow=numdraws, ncol=J)
for(j in 1:J) { # scale by tau * lambda
    beta.sim[, j] <- beta.sim[, j] * tau * lambda_tilde[, j]
}

## simulate observed gamma_hat estimates conditional on beta.sim
gamma_hat <- matrix(rnorm(n=numdraws * J), nrow=numdraws, ncol=J)
for(j in 1:J) { # scale by se.gamma_hat and shift by beta.sim
    gamma_hat[, j] <- beta.sim[, j] + gamma_hat[, j] * se.gamma_hat[j]
}

## simulate observed alpha_hat estimates conditional on alpha.sim
alpha_hat <- matrix(rnorm(n=numdraws * J), nrow=numdraws, ncol=J)
for(j in 1:J) { # scale by se.alpha_hat and shift by alpha.sim
    alpha_hat[, j] <- alpha.sim[, j] + alpha_hat[, j] * se.alpha_hat[j]
}

## simulate coefficient ratios, possibly estimated by delta method
thetaIV.sim <- gamma_hat / alpha_hat + ifelse(use.delta,
                                              se.alpha_hat^2 * gamma_hat / alpha_hat^3,
                                              0)

se.thetaIV.sim <- sqrt((se.gamma_hat / alpha_hat)^2 +
                       gamma_hat^2 * se.alpha_hat^2 / alpha_hat^4)

inv.var <- 1 / se.thetaIV.sim^2

theta_IVW <- rowSums(thetaIV.sim * inv.var) / rowSums(inv.var) 
se.theta_IVW  <- sqrt(1 / rowSums(inv.var))
z.theta_IVW <- theta_IVW / se.theta_IVW

## calculate


wmed <- numeric(numdraws)
wmed.sd.bowden <- numeric(numdraws)

for(i in 1:numdraws) {
    wmed[i] <- matrixStats::weightedMedian(thetaIV.sim[i, ], w=inv.var[i, ])

    wmed.sd.bowden[i] <- weighted.median.sd.bowden(alpha_hat[i, ], gamma_hat[i, ],
                                                   se.alpha_hat,
                                                   se.gamma_hat, weights=inv.var[i, ])
}

z.theta_IVW.bowden <- wmed / wmed.sd.bowden

sd(z.theta_IVW.bowden)
