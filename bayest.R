# bayest - Bayesian t-tests
#     Copyright (C) 2020  Riko Kelter
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

bayes.t.test <-
  function(n = 10000,
           plot = "all",
           firstComp,
           secondComp,
           hyperpars = "wide",
           ci = "0.95",
           burnin = n / 2,
           sd = "sd",
           q = 0.1) {
    # utility function for computing the mode
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    # Set the credible level
    credibleLevel = ci
    
    Nsim = n # number of steps performed by the Gibbs sampler
    firstComponent = firstComp # set first component
    secondComponent = secondComp # set second component
    mu1 <- numeric(Nsim) # arrays for parameters in first chain
    mu2 <- numeric(Nsim)
    sigma1Sq <- numeric(Nsim)
    sigma2Sq <- numeric(Nsim)
    
    mu1SecC <- numeric(Nsim) # arrays for parameters in first chain
    mu2SecC <- numeric(Nsim)
    sigma1SqSecC <- numeric(Nsim)
    sigma2SqSecC <- numeric(Nsim)
    
    mu1ThdC <- numeric(Nsim) # arrays for parameters in first chain
    mu2ThdC <- numeric(Nsim)
    sigma1SqThdC <- numeric(Nsim)
    sigma2SqThdC <- numeric(Nsim)
    
    mu1FrtC <- numeric(Nsim) # arrays for parameters in first chain
    mu2FrtC <- numeric(Nsim)
    sigma1SqFrtC <- numeric(Nsim)
    sigma2SqFrtC <- numeric(Nsim)
    
    mu1[1] = mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
    mu2[1] = mean(secondComponent)
    sigma1Sq[1] = var(firstComponent)
    sigma2Sq[1] = var(secondComponent)
    
    mu1SecC[1] = 10 * mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
    mu2SecC[1] = 10 * mean(secondComponent)
    sigma1SqSecC[1] = 10 * var(firstComponent)
    sigma2SqSecC[1] = 10 * var(secondComponent)
    
    mu1ThdC[1] = 1 / 10 * mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
    mu2ThdC[1] = 1 / 10 * mean(secondComponent)
    sigma1SqThdC[1] = 1 / 10 * var(firstComponent)
    sigma2SqThdC[1] = 1 / 10 * var(secondComponent)
    
    mu1FrtC[1] = 5 * mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
    mu2FrtC[1] = 5 * mean(secondComponent)
    sigma1SqFrtC[1] = 5 * var(firstComponent)
    sigma2SqFrtC[1] = 5 * var(secondComponent)
    
    # Formula parts
    N_1_S = length(firstComponent)
    N_2_S = length(secondComponent)
    vary2 = var(secondComponent)
    vary1 = var(firstComponent)
    bary1 = mean(firstComponent)
    bary2 = mean(secondComponent)
    
    # Set hyperparameters
    if (hyperpars == "rafterys") {
      smple = c(firstComponent, secondComponent)
      b0 = mean(smple) # Raftery's hyperparameters
      B0 = var(smple) / 2.6 * (max(smple) - min(smple))
      #B0=sd(smple)*10
      c0 = 1.28 # Raftery's hyperparameters
      C0 = 0.36 * var(smple) # Raftery's hyperparameters
      #C0=10*sd(smple)
    }
    if (hyperpars == "custom") {
      smple = c(firstComponent, secondComponent)
      b0 = mean(smple) # Raftery's hyperparameters
      B0 = 250 * var(smple)
      c0 = 1.28 # Raftery's hyperparameters
      C0 = -0.5 * var(smple) + (c0 + 0.5 * (length(firstComp) + length(secondComp))) *
        q
      
    }
    if (hyperpars == "wide") {
      smple = c(firstComponent, secondComponent)
      b0 = mean(smple)
      B0 = 10 * var(smple)
      c0 = 0.01
      C0 = 0.01
    }
    if (hyperpars == "medium") {
      smple = c(firstComponent, secondComponent)
      b0 = mean(smple)
      B0 = 5 * var(smple)
      c0 = 0.1
      C0 = 0.1
    }
    if (hyperpars == "narrow") {
      smple = c(firstComponent, secondComponent)
      b0 = mean(smple)
      B0 = 1 * var(smple)
      c0 = 1
      C0 = 1
    }
    
    
    # Gibbs sampling via full conditionals
    #library(MCMCpack)
    for (t in 2:Nsim) {
      # FIRST CHAIN
      # sample sigma_1^2|mu1,mu2,sigma_2^2,S,y
      c_1_S = c0 + 0.5 * N_1_S
      C_1_S = C0 + 0.5 * sum((firstComponent - mu1[t - 1]) ^ 2)
      sigma1Sq[t] = MCMCpack::rinvgamma(1, shape = c_1_S, scale = C_1_S)
      
      # sample sigma_2^2|mu1,mu2,sigma_1^2,S,y
      c_2_S = c0 + 0.5 * N_2_S
      C_2_S = C0 + 0.5 * sum((secondComponent - mu2[t - 1]) ^ 2)
      sigma2Sq[t] = MCMCpack::rinvgamma(1, shape = c_2_S, scale = C_2_S)
      
      # sample mu1|mu2,sigma_1^2,sigma_2^2,S,y
      B_1_S = 1 / ((1 / B0) + (1 / sigma1Sq[t]) * N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
      b_1_S = B_1_S * ((1 / sigma1Sq[t]) * N_1_S * bary1 + (1 / B0) * b0)
      mu1[t] = rnorm(1, mean = b_1_S, sd = B_1_S)
      
      # sample mu2|mu1,sigma_1^2,sigma_2^2,S,y
      B_2_S = 1 / ((1 / B0) + (1 / sigma2Sq[t]) * N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
      b_2_S = B_2_S * ((1 / sigma2Sq[t]) * N_2_S * bary2 + (1 / B0) * b0)
      mu2[t] = rnorm(1, mean = b_2_S, sd = B_2_S)
      
      ###############################################################
      # SECOND CHAIN
      # sample sigma_1^2|mu1,mu2,sigma_2^2,S,y
      C_1_S = C0 + 0.5 * sum((firstComponent - mu1SecC[t - 1]) ^ 2)
      sigma1SqSecC[t] = MCMCpack::rinvgamma(1, shape = c_1_S, scale = C_1_S)
      
      # sample sigma_2^2|mu1,mu2,sigma_1^2,S,y
      C_2_S = C0 + 0.5 * sum((secondComponent - mu2SecC[t - 1]) ^ 2)
      sigma2SqSecC[t] = MCMCpack::rinvgamma(1, shape = c_2_S, scale = C_2_S)
      
      # sample mu1|mu2,sigma_1^2,sigma_2^2,S,y
      B_1_S = 1 / ((1 / B0) + (1 / sigma1SqSecC[t]) * N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
      b_1_S = B_1_S * ((1 / sigma1SqSecC[t]) * N_1_S * bary1 + (1 / B0) *
                         b0)
      mu1SecC[t] = rnorm(1, mean = b_1_S, sd = B_1_S)
      
      # sample mu2|mu1,sigma_1^2,sigma_2^2,S,y
      B_2_S = 1 / ((1 / B0) + (1 / sigma2SqSecC[t]) * N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
      b_2_S = B_2_S * ((1 / sigma2SqSecC[t]) * N_2_S * bary2 + (1 / B0) *
                         b0)
      mu2SecC[t] = rnorm(1, mean = b_2_S, sd = B_2_S)
      
      ###############################################################
      # THIRD CHAIN
      # sample sigma_1^2|mu1,mu2,sigma_2^2,S,y
      C_1_S = C0 + 0.5 * sum((firstComponent - mu1ThdC[t - 1]) ^ 2)
      sigma1SqThdC[t] = MCMCpack::rinvgamma(1, shape = c_1_S, scale = C_1_S)
      
      # sample sigma_2^2|mu1,mu2,sigma_1^2,S,y
      C_2_S = C0 + 0.5 * sum((secondComponent - mu2ThdC[t - 1]) ^ 2)
      sigma2SqThdC[t] = MCMCpack::rinvgamma(1, shape = c_2_S, scale = C_2_S)
      
      # sample mu1|mu2,sigma_1^2,sigma_2^2,S,y
      B_1_S = 1 / ((1 / B0) + (1 / sigma1SqThdC[t]) * N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
      b_1_S = B_1_S * ((1 / sigma1SqThdC[t]) * N_1_S * bary1 + (1 / B0) *
                         b0)
      mu1ThdC[t] = rnorm(1, mean = b_1_S, sd = B_1_S)
      
      # sample mu2|mu1,sigma_1^2,sigma_2^2,S,y
      B_2_S = 1 / ((1 / B0) + (1 / sigma2SqThdC[t]) * N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
      b_2_S = B_2_S * ((1 / sigma2SqThdC[t]) * N_2_S * bary2 + (1 / B0) *
                         b0)
      mu2ThdC[t] = rnorm(1, mean = b_2_S, sd = B_2_S)
      
      ###############################################################
      # FOURTH CHAIN
      # sample sigma_1^2|mu1,mu2,sigma_2^2,S,y
      C_1_S = C0 + 0.5 * sum((firstComponent - mu1FrtC[t - 1]) ^ 2)
      sigma1SqFrtC[t] = MCMCpack::rinvgamma(1, shape = c_1_S, scale = C_1_S)
      
      # sample sigma_2^2|mu1,mu2,sigma_1^2,S,y
      C_2_S = C0 + 0.5 * sum((secondComponent - mu2FrtC[t - 1]) ^ 2)
      sigma2SqFrtC[t] = MCMCpack::rinvgamma(1, shape = c_2_S, scale = C_2_S)
      
      # sample mu1|mu2,sigma_1^2,sigma_2^2,S,y
      B_1_S = 1 / ((1 / B0) + (1 / sigma1SqFrtC[t]) * N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
      b_1_S = B_1_S * ((1 / sigma1SqFrtC[t]) * N_1_S * bary1 + (1 / B0) *
                         b0)
      mu1FrtC[t] = rnorm(1, mean = b_1_S, sd = B_1_S)
      
      # sample mu2|mu1,sigma_1^2,sigma_2^2,S,y
      B_2_S = 1 / ((1 / B0) + (1 / sigma2SqFrtC[t]) * N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
      b_2_S = B_2_S * ((1 / sigma2SqFrtC[t]) * N_2_S * bary2 + (1 / B0) *
                         b0)
      mu2FrtC[t] = rnorm(1, mean = b_2_S, sd = B_2_S)
      
      ##############################################################
      # All chains samples, reiterate
    }
    
    # Adapt chains for selected burnin
    mu1 = mu1[burnin:Nsim]
    mu1SecC = mu1SecC[burnin:Nsim]
    mu1ThdC = mu1ThdC[burnin:Nsim]
    mu1FrtC = mu1FrtC[burnin:Nsim]
    
    mu2 = mu2[burnin:Nsim]
    mu2SecC = mu2SecC[burnin:Nsim]
    mu2ThdC = mu2ThdC[burnin:Nsim]
    mu2FrtC = mu2FrtC[burnin:Nsim]
    
    sigma1Sq = sigma1Sq[burnin:Nsim]
    sigma1SqSecC = sigma1SqSecC[burnin:Nsim]
    sigma1SqThdC = sigma1SqThdC[burnin:Nsim]
    sigma1SqFrtC = sigma1SqFrtC[burnin:Nsim]
    
    sigma2Sq = sigma2Sq[burnin:Nsim]
    sigma2SqSecC = sigma2SqSecC[burnin:Nsim]
    sigma2SqThdC = sigma2SqThdC[burnin:Nsim]
    sigma2SqFrtC = sigma2SqFrtC[burnin:Nsim]
    
    if (sd == "sd") {
      sigma1Sq = sqrt(sigma1Sq)
      sigma1SqSecC = sqrt(sigma1SqSecC)
      sigma1SqThdC = sqrt(sigma1SqThdC)
      sigma1SqFrtC = sqrt(sigma1SqFrtC)
      
      sigma2Sq = sqrt(sigma2Sq)
      sigma2SqSecC = sqrt(sigma2SqSecC)
      sigma2SqThdC = sqrt(sigma2SqThdC)
      sigma2SqFrtC = sqrt(sigma2SqFrtC)
    }
    
    diffOfMeans = mu2 - mu1
    diffOfMeansSecC = mu2SecC - mu1SecC
    diffOfMeansThdC = mu2ThdC - mu1ThdC
    diffOfMeansFrtC = mu2FrtC - mu1FrtC
    
    diffOfVariances = sigma2Sq - sigma1Sq
    diffOfVariancesSecC = sigma2SqSecC - sigma1SqSecC
    diffOfVariancesThdC = sigma2SqThdC - sigma1SqThdC
    diffOfVariancesFrtC = sigma2SqFrtC - sigma1SqFrtC
    
    # Effect size
    if (sd == "var") {
      effectSize = (mu2 - mu1) / (sqrt((sigma1Sq + sigma2Sq) / 2))
      effectSizeSecC = (mu2SecC - mu1SecC) / (sqrt((sigma1SqSecC + sigma2SqSecC) /
                                                     2))
      effectSizeThdC = (mu2ThdC - mu1ThdC) / (sqrt((sigma1SqThdC + sigma2SqThdC) /
                                                     2))
      effectSizeFrtC = (mu2FrtC - mu1FrtC) / (sqrt((sigma1SqFrtC + sigma2SqFrtC) /
                                                     2))
    }
    if (sd == "sd") {
      effectSize = (mu2 - mu1) / (sqrt((sigma1Sq ^ 2 + sigma2Sq ^ 2) / 2))
      effectSizeSecC = (mu2SecC - mu1SecC) / (sqrt((sigma1SqSecC ^ 2 + sigma2SqSecC ^
                                                      2) / 2))
      effectSizeThdC = (mu2ThdC - mu1ThdC) / (sqrt((sigma1SqThdC ^ 2 + sigma2SqThdC ^
                                                      2) / 2))
      effectSizeFrtC = (mu2FrtC - mu1FrtC) / (sqrt((sigma1SqFrtC ^ 2 + sigma2SqFrtC ^
                                                      2) / 2))
    }
    
    # output
    if (plot == "all") {
      plotpar <- par(mfrow = c(2, 2))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(
        mu1,
        freq = FALSE,
        main = "",
        col = "cornflowerblue",
        border = "white",
        xlab = expression(mu[1])
      )
      plot(
        mu1,
        ty = "l",
        col = "cornflowerblue",
        xlab = expression(mu[1]),
        ylab = ""
      )
      mcmcMatrixMu1 = matrix(0,
                             nrow = length(mu1),
                             ncol = 4,
                             byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = mu1)
      mcmcObj2 = coda::mcmc(data = mu1SecC)
      mcmcObj3 = coda::mcmc(data = mu1ThdC)
      mcmcObj4 = coda::mcmc(data = mu1FrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(mu1)
      
      
      
      plotpar <- par(mfrow = c(2, 2))
      hist(
        mu2,
        freq = FALSE,
        main = "",
        col = "cornflowerblue",
        border = "white",
        xlab = expression(mu[2])
      )
      plot(
        mu2,
        ty = "l",
        col = "cornflowerblue",
        xlab = expression(mu[2]),
        ylab = ""
      )
      mcmcMatrixMu2 = matrix(0,
                             nrow = length(mu2),
                             ncol = 4,
                             byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = mu2)
      mcmcObj2 = coda::mcmc(data = mu2SecC)
      mcmcObj3 = coda::mcmc(data = mu2ThdC)
      mcmcObj4 = coda::mcmc(data = mu2FrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(mu2)
      
      plotpar <- par(mfrow = c(2, 2))
      if (sd == "var") {
        hist(
          sigma1Sq,
          freq = FALSE,
          main = "",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[1] ^ 2)
        )
        plot(
          sigma1Sq,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[1] ^ 2),
          ylab = ""
        )
      }
      if (sd == "sd") {
        hist(
          sigma1Sq,
          freq = FALSE,
          main = "",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[1])
        )
        plot(
          sigma1Sq,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[1]),
          ylab = ""
        )
      }
      mcmcMatrixSigma1Sq = matrix(0,
                                  nrow = length(sigma1Sq),
                                  ncol = 4,
                                  byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = sigma1Sq)
      mcmcObj2 = coda::mcmc(data = sigma1SqSecC)
      mcmcObj3 = coda::mcmc(data = sigma1SqThdC)
      mcmcObj4 = coda::mcmc(data = sigma1SqFrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(sigma1Sq)
      
      
      plotpar <- par(mfrow = c(2, 2))
      if (sd == "var") {
        hist(
          sigma2Sq,
          freq = FALSE,
          main = "",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[2] ^ 2)
        )
        plot(
          sigma2Sq,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[2] ^ 2),
          ylab = ""
        )
      }
      if (sd == "sd") {
        hist(
          sigma2Sq,
          freq = FALSE,
          main = "",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[2])
        )
        plot(
          sigma2Sq,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[2]),
          ylab = ""
        )
      }
      mcmcMatrixSigma2Sq = matrix(0,
                                  nrow = length(sigma2Sq),
                                  ncol = 4,
                                  byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = sigma2Sq)
      mcmcObj2 = coda::mcmc(data = sigma2SqSecC)
      mcmcObj3 = coda::mcmc(data = sigma2SqThdC)
      mcmcObj4 = coda::mcmc(data = sigma2SqFrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(sigma2Sq)
      
      plotpar <- par(mfrow = c(2, 2))
      hist(
        diffOfMeans,
        freq = FALSE,
        main = "Difference of means",
        col = "cornflowerblue",
        border = "white",
        xlab = expression(mu[2] - mu[1])
      )
      plot(
        diffOfMeans,
        ty = "l",
        col = "cornflowerblue",
        xlab = expression(mu[2] - mu[1]),
        ylab = ""
      )
      mcmcMatrixDiffOfMeans = matrix(0,
                                     nrow = length(diffOfMeans),
                                     ncol = 4,
                                     byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = diffOfMeans)
      mcmcObj2 = coda::mcmc(data = diffOfMeansSecC)
      mcmcObj3 = coda::mcmc(data = diffOfMeansThdC)
      mcmcObj4 = coda::mcmc(data = diffOfMeansFrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeans)
      
      plotpar <- par(mfrow = c(2, 2))
      if (sd == "var") {
        hist(
          diffOfVariances,
          freq = FALSE,
          main = "Difference of variances",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[2] ^ 2 - sigma[1] ^ 2)
        )
        plot(
          diffOfVariances,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[2] ^ 2 - sigma[1] ^ 2),
          ylab = ""
        )
      }
      if (sd == "sd") {
        hist(
          diffOfVariances,
          freq = FALSE,
          main = "Difference of standard deviations",
          col = "cornflowerblue",
          border = "white",
          xlab = expression(sigma[2] - sigma[1])
        )
        plot(
          diffOfVariances,
          ty = "l",
          col = "cornflowerblue",
          xlab = expression(sigma[2] - sigma[1]),
          ylab = ""
        )
      }
      mcmcMatrixDiffOfVariances = matrix(
        0,
        nrow = length(diffOfVariances),
        ncol = 4,
        byrow = FALSE
      )
      mcmcObj1 = coda::mcmc(data = diffOfVariances)
      mcmcObj2 = coda::mcmc(data = diffOfVariancesSecC)
      mcmcObj3 = coda::mcmc(data = diffOfVariancesThdC)
      mcmcObj4 = coda::mcmc(data = diffOfVariancesFrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariances)
      
      plotpar <- par(mfrow = c(2, 2))
      hist(
        effectSize,
        freq = FALSE,
        main = "Effect size",
        col = "cornflowerblue",
        border = "white",
        xlab = expression(delta)
      )
      plot(
        effectSize,
        ty = "l",
        col = "cornflowerblue",
        xlab = expression(delta),
        ylab = ""
      )
      mcmcMatrixEffectSize = matrix(0,
                                    nrow = length(effectSize),
                                    ncol = 4,
                                    byrow = FALSE)
      mcmcObj1 = coda::mcmc(data = effectSize)
      mcmcObj2 = coda::mcmc(data = effectSizeSecC)
      mcmcObj3 = coda::mcmc(data = effectSizeThdC)
      mcmcObj4 = coda::mcmc(data = effectSizeFrtC)
      mcmcListObj = coda::mcmc.list(mcmcObj1, mcmcObj2, mcmcObj3, mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize)
      
      # Set 2x1 par
      oldpar <- par(mfrow = c(2, 1))
      
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(
        effectSize,
        freq = FALSE,
        col = "cornflowerblue",
        border = "white",
        xlab = expression(delta),
        main = mtext(
          "Posterior distribution of Effect size " ~ delta,
          side = 3,
          line = 2.2,
          col = "black",
          cex = 1,
          font = 2
        ),
        xaxt = "n"
      )
      ticks <- seq(from = -3, to = 3, by = 0.1)
      # And draw the axis:
      axis(1, at = ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(
        x = c(getmode(effectSize), getmode(effectSize)),
        y = c(-0.025, 0.025),
        type = "l",
        lwd = 2.5,
        col = "blue"
      )
      
      # Posterior CI
      postCI = quantile(effectSize, probs = c((1 - ci) / 2, ci + (1 - ci) /
                                                2))
      lines(
        x = c(postCI[1], postCI[2]),
        y = c(0, 0),
        type = "l",
        lwd = 2,
        col = "black"
      )
      lines(
        x = c(postCI[1], postCI[1]),
        y = c(-0.02, 0.02),
        type = "l",
        lwd = 2,
        col = "black"
      )
      lines(
        x = c(postCI[2], postCI[2]),
        y = c(-0.02, 0.02),
        type = "l",
        lwd = 2,
        col = "black"
      )
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(
        x = c(-0.2, 0.2),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "red"
      )
      # ROPE for small effect
      lines(
        x = c(-0.5,-0.2),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "orange"
      )
      lines(
        x = c(0.2, 0.5),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "orange"
      )
      
      # ROPE for medium effect
      lines(
        x = c(-0.8,-0.5),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "green"
      )
      lines(
        x = c(0.5, 0.8),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "green"
      )
      # ROPE for large effect
      lines(
        x = c(-2,-0.8),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "purple"
      )
      lines(
        x = c(0.8, 2),
        y = c(0.04, 0.04),
        type = "l",
        lwd = 2,
        col = "purple"
      )
      
      abline(
        v = 0.8,
        col = "purple",
        lwd = 0.5,
        lty = "dashed"
      )
      abline(
        v = -0.8,
        col = "purple",
        lwd = 0.5,
        lty = "dashed"
      )
      abline(
        v = 0.5,
        col = "green",
        lwd = 0.5,
        lty = "dashed"
      )
      abline(
        v = -0.5,
        col = "green",
        lwd = 0.5,
        lty = "dashed"
      )
      abline(
        v = 0.2,
        col = "orange",
        lwd = 0.5,
        lty = "dashed"
      )
      abline(
        v = -0.2,
        col = "orange",
        lwd = 0.5,
        lty = "dashed"
      )
      
      # Text above the plot
      mtext(
        paste0(
          ci * 100,
          "% ",
          "CI: [",
          round(postCI[1], 3),
          ",",
          round(postCI[2], 3),
          "]"
        ),
        side = 3,
        line = 0.15,
        at = min(effectSize),
        col = "black",
        cex = 1,
        font = 2
      )
      mtext(
        paste0("Post. Mean: ", round(mean(effectSize), 3)),
        side = 3,
        line = 1.15,
        at = max(effectSize),
        col = "black",
        cex = 1,
        font = 2
      )
      mtext(
        paste0("Post. Mode: ", round(getmode(effectSize), 3)),
        side = 3,
        line = 0.15,
        at = max(effectSize),
        col = "black",
        cex = 1,
        font = 2
      )
      
      
      # Bar plot
      effSizeCI = quantile(effectSize, probs = c(((1 - ci) / 2), (ci + (1 -
                                                                          ci) / 2)))
      effectSizeCIValues = effectSize[effectSize > effSizeCI[1] &
                                        effectSize < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8)) /
        length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8)) /
        length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 &
                                                 effectSizeCIValues < 0.8)) / length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 &
                                                 effectSizeCIValues > -0.8)) / length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 &
                                                effectSizeCIValues < 0.5)) / length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 &
                                                effectSizeCIValues > -0.5)) / length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 &
                                          effectSizeCIValues > -0.2)) / length(effectSizeCIValues)
      
      lbls = c(
        "large positive:",
        "medium positive:",
        "small positive:",
        "no effect:",
        "small negative:",
        "medium negative:",
        "large negative:"
      )
      pct = c(
        largeNegEffectIterations,
        mediumNegEffectIterations,
        smallNegEffectIterations,
        noEffectIterations,
        smallPosEffectIterations,
        mediumPosEffectIterations,
        largePosEffectIterations
      )
      pct = round(pct, digits = 4)
      pct = pct * 100
      lbls = paste(lbls, pct)
      lbls <- paste(pct, "%", sep = "") # ad % to labels
      
      # Barplot
      b <- barplot(
        pct,
        main = "Posterior ROPE-Analysis",
        xlab = "ROPEs of standard effect sizes",
        ylab = "Posterior-Percentage included inside ROPE",
        names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
        col = "cornflowerblue",
        border = "white",
        horiz = FALSE,
        ylim = c(0, 107.5)
      )
      b
      text(x = b,
           y = pct + 4,
           labels = lbls)
      
      # print out modes
      getmode(diffOfMeans)
      getmode(diffOfVariances)
      
      # return object of both modes
      posteriorModes = c(getmode(diffOfMeans), getmode(diffOfVariances))
      posteriorModes
    }
    if (plot == "none") {
      # return dataframe of both modes, do not print/plot anything
      if (sd == "var") {
        Parameter = c(
          "Difference of means (mu2-mu1)",
          "Difference of variances (sigma2^2-sigma1^2)",
          "Effect size",
          "MAPE"
        )
      }
      if (sd == "sd") {
        Parameter = c(
          "Difference of means (mu2-mu1)",
          "Difference of standard deviations (sigma2-sigma1)",
          "Effect size",
          "MAPE"
        )
      }
      PosteriorMode = c(getmode(diffOfMeans),
                        getmode(diffOfVariances),
                        getmode(effectSize),-1)
      PosteriorExpectation = c(mean(diffOfMeans),
                               mean(diffOfVariances),
                               mean(effectSize),-1)
      LowerCI = c(
        quantile(diffOfMeans, probs = c((1 - credibleLevel) / 2, credibleLevel +
                                          (1 - credibleLevel) / 2
        ))[1],
        quantile(diffOfVariances, probs = c((1 - credibleLevel) / 2, credibleLevel +
                                              (1 - credibleLevel) / 2
        ))[1],
        quantile(effectSize, probs = c((1 - credibleLevel) / 2, credibleLevel +
                                         (1 - credibleLevel) / 2
        ))[1],-1
      )
      UpperCI = c(
        quantile(diffOfMeans, probs = c(1 - credibleLevel, credibleLevel))[2],
        quantile(diffOfVariances, probs = c(1 - credibleLevel, credibleLevel))[2],
        quantile(effectSize, probs = c((1 - credibleLevel) / 2, credibleLevel +
                                         (1 - credibleLevel) / 2
        ))[2],-1
      )
      
      effSizeCI = quantile(effectSize, probs = c(((1 - ci) / 2), (ci + (1 -
                                                                          ci) / 2)))
      effectSizeCIValues = effectSize[effectSize > effSizeCI[1] &
                                        effectSize < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8)) /
        length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8)) /
        length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 &
                                                 effectSizeCIValues < 0.8)) / length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 &
                                                 effectSizeCIValues > -0.8)) / length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 &
                                                effectSizeCIValues < 0.5)) / length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 &
                                                effectSizeCIValues > -0.5)) / length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 &
                                          effectSizeCIValues > -0.2)) / length(effectSizeCIValues)
      
      pct = c(
        largeNegEffectIterations,
        mediumNegEffectIterations,
        smallNegEffectIterations,
        noEffectIterations,
        smallPosEffectIterations,
        mediumPosEffectIterations,
        largePosEffectIterations
      )
      pct = round(pct, digits = 4)
      pct = pct * 100
      mx = max(pct)
      indx = which(pct == mx)
      if (indx == 1) {
        str <- "Large Negative"
      }
      if (indx == 2) {
        str <- "Medium Negative"
      }
      if (indx == 3) {
        str <- "Small Negative"
      }
      if (indx == 4) {
        str <- "No Effect"
      }
      if (indx == 5) {
        str <- "Small Positive"
      }
      if (indx == 6) {
        str <- "Medium Positive"
      }
      if (indx == 7) {
        str <- "Large Positive"
      }
      MAPE <- c(0, 0, 0, indx)
      df = data.frame(Parameter,
                      PosteriorMode,
                      PosteriorExpectation,
                      LowerCI,
                      UpperCI,
                      MAPE)
      df
    }
  }
