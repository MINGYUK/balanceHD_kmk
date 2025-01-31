rm(list = ls())

library(balanceHDkmk)
library(ebalTorchkmk)
library(mvtnorm)
library(Matrix)

# Read inputs to simulation

args=(commandArgs(TRUE))
beta.setup = 1
prop.setup = 1
n = 500
p = 10
eps = 0.2
C = 1
extra.param = 1

NREP = 1

experiment = 1

tau = 1

run.comparison = function(XYW) {
  
  X = XYW$X
  Y = XYW$Y
  W = XYW$W
  
  # Run residual balancing
  tau.rb = residualBalance.ate(X, Y, W, target.pop = c(0, 1), fit.method = "elnet", alpha = 0.9, zeta = 0.5)
  tau.rb.plain = residualBalance.ate(X, Y, W, target.pop = c(0, 1), fit.method = "none", zeta = 0.5)
  
  # Run inverse-propensity weighted methods
  tau.ipw.elnet = ipw.ate(X, Y, W, target.pop = 1, fit.method = "none", prop.method = "elnet", alpha.prop = 0.5)
  tau.aipw = ipw.ate(X, Y, W, target.pop = 1, prop.weighted.fit = FALSE, targeting.method = "AIPW",
                     fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)
  tau.ipw.weighted = ipw.ate(X, Y, W, target.pop = 1, prop.weighted.fit = TRUE, targeting.method = "AIPW",
                             fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)
  tau.tmle = ipw.ate(X, Y, W, prop.weighted.fit = FALSE, targeting.method = "TMLE",
                     target.pop = 1, fit.method = "elnet", prop.method = "elnet", alpha.fit = 0.9, alpha.prop = 0.5)
  
  ################### Custom implementation by manjimin
  # Run EBAL
  tau.ebal = ebal.ate(Treatment = W,
                      X = X,
                      Y=Y,
                      method='AutoDiff',
                      constraint.tolerance = 1)
  
  # Run SBW
  tau.sbw = sbw.ate(W,
                    X,
                    Y,
                    tolerance = 0.02,
                    algorithm = FALSE,
                    std = 'group',
                    solver = 'osqp')
  
  # Run dynbal
  tau.dynbal = dynbal.ate(Treatment = W,
                          X = X,
                          Y=Y,
                          method='AutoDiff',
                          constraint.tolerance = 1)
  
  #####################################################
  # Run other baselines
  tau.naive = naive.ate(X, Y, W)
  tau.elnet = elnet.ate(X, Y, W, target.pop = 1, alpha = 0.9)
  tau.twostep = twostep.lasso.ate(X, Y, W, target.pop = 1)
  
  results = c(Naive=tau.naive,
              Elnet=tau.elnet,
              Twostep = tau.twostep,
              ResidBalance=tau.rb,
              ResidBalancePlain=tau.rb.plain,
              IPW.Elnet = tau.ipw.elnet,
              IPW.Residual = tau.aipw,
              IPW.Weighted = tau.ipw.weighted,
              TMLE = tau.tmle,
              EBAL = tau.ebal,
              SBW = tau.sbw,
              DYNBAL = tau.dynbal
  )
  return(results)
  
}

# Intialize beta

if(beta.setup == 0) {
  beta.raw = rep(1, p)
} else if(beta.setup == 1){
  beta.raw = 1/sqrt(1:p)
} else if(beta.setup == 2) {
  beta.raw = 1/(9 + 1:p)
} else if(beta.setup == 3) {
  beta.raw = c(rep(10, 10), rep(1, 90), rep(0, p - 100))
} else if(beta.setup == 4) {
  beta.raw = c(rep(10, 10), rep(0, p - 10))
} else if(beta.setup == 5) {
  beta.raw = 1/(1:p)^2
} else if(beta.setup == 6) {
  beta.raw = 1/(1:p)
}

beta.main =  C  * beta.raw / sqrt(sum(beta.raw^2))

# Initialize data-generating model
if(prop.setup == 1) {
  
  # simple setup
  if (extra.param == 1) {
    delta.clust = 4 / sqrt(n) * rep(1, p)
  } else if (extra.param == 2) {
    delta.clust = 40 / sqrt(n) * rep(c(1, rep(0, 9)), p/10)
  }
  gen.data = function() {
    CLUST = rbinom(n, 1, 0.5)
    W = rbinom(n, 1, eps + (1 - 2*eps) * CLUST)
    X = matrix(rnorm(n * p), n, p) + outer(CLUST, delta.clust)
    Y = X %*% beta.main + rnorm(n) + tau * W
    list(X=X, Y=Y, W=W, catt=tau)
  }
  
} else if(prop.setup == 2) {
  
  # AR-1 model
  rho = c(0.5, 0.9)[extra.param]
  beta.main = beta.main[1 + (23 * (0:(length(beta.main) - 1))) %% length(beta.main)]
  beta.prop = c(rep(1/40, 100), rep(0, p - 100))
  gen.data = function() {
    X = 10 * rmvnorm(n, sigma = outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y))), method = "chol")
    W = rbinom(n, 1, prob = 1/(1 + exp(-beta.prop * X)))
    Y = X %*% beta.main + rnorm(n, 0, 1) + tau * W
    list(X=X, Y=Y, W=W, catt=tau)
  }
  
} else if(prop.setup == 3) {
  
  # many cluster model
  nclust = c(6, 20)[extra.param]
  beta.main = sqrt(nclust) * beta.main
  clust.ptreat = rep(c(eps, 1 - eps), nclust/2)
  tau.clust = tau * rexp(nclust)
  gen.data = function() {
    cluster.center = matrix(rnorm(nclust * p), nclust, p)
    cluster = sample.int(nclust, n, replace = TRUE)
    X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
    W = rbinom(n, 1, clust.ptreat[cluster])
    Y = X %*% beta.main + rnorm(n, 0, 1) + tau.clust[cluster] * W
    list(X=X, Y=Y, W=W, catt=mean(tau.clust[cluster[W==1]]))
  }
  
} else if(prop.setup == 4) {
  
  tau = 1
  
  # misspecified model
  gen.data = function() {
    X = matrix(rnorm(n * p), n, p)
    # Normalized such that E[tauX | W = 1] = 1
    tauX = log(1 + exp(-2 - 2 * X[,1])) / 0.915
    ptreat = 1 - exp(-tauX)
    W = rbinom(n, 1, prob = ptreat)
    mean(tauX[W == 1])
    Y = rnorm(n, 0, 1) + tauX * (2 * W - 1) / 2 + rowSums(X[,1:10])
    list(X=X, Y=Y, W=W, catt=mean(tauX[W==1]))
  }
  
} else if(prop.setup == 6) {
  
  tau = 0.5
  rho = 0.5
  
  Sigma = outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  
  beta.raw = 1/(1:p)^2
  beta.main =  C  * beta.raw / sqrt(sum(beta.raw^2))
  
  if (extra.param == 1) {
    beta.prop = 1/(1:p)^2
  } else if (extra.param == 2) {
    beta.prop = 1/sqrt(1:p)
  }
  
  beta.prop = eps * beta.prop / sqrt(sum(beta.prop^2))
  
  gen.data = function() {
    X = rmvnorm(n, sigma = Sigma, method = "chol")
    theta = X %*% beta.prop + rnorm(n)
    W = rbinom(n, 1, 1/(1 + exp(theta)))
    Y = X %*% beta.main + W * tau + rnorm(n)
    list(X=X, Y=Y, W=W, catt=tau)
  }
  
}

gen.data = function() {
  CLUST = rbinom(n, 1, 0.5)
  delta.clust = 4 / sqrt(n) * rep(1, p)
  W = rbinom(n, 1, eps + (1 - 2*eps) * CLUST)
  X = matrix(rnorm(n * p), n, p) + outer(CLUST, delta.clust)
  Y = X %*% beta.main + rnorm(n) + tau * W
  list(X=X, Y=Y, W=W, catt=tau)
}

data = gen.data()
X = data$X
Y = data$Y
W = data$W


if (experiment == 1) {
  
  # MSE vs. other methods
  results = replicate(NREP, {
    XYW = gen.data()
    res = run.comparison(XYW)
    print(c(CATT=XYW$catt, res))
    return(c(CATT=XYW$catt, res))
  })
  
  print(colMeans((t(results) - results["CATT",])/results["CATT",]))
  print(colMeans(((t(results) - results["CATT",])/results["CATT",])^2))
  
} else if (experiment == 2) {
  
  # coverage
  results = replicate(NREP, {
    XYW = gen.data()
    res = coverage.test(XYW)
    print(c(CATT=XYW$catt, res))
    return(c(CATT=XYW$catt, res))
  })
  rownames(results) = c("CATT", "tau.hat", "se.hat")
  se.err = (results[2,] - results[1,]) / results[3,]
  coverage = mean(abs(se.err) <= qnorm(0.975))
}

library(ebalTorchkmk)
library(Matrix)


save.image(file=paste0(paste("results/res", beta.setup, prop.setup, n, p, eps, C, extra.param, NREP, experiment, sep = "-"), ".RData"))
