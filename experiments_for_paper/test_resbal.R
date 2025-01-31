library(Matrix)
library(balanceHD)
# Hyperparameters
n = 500
p = 50000
eps = 0.5
C = 1

tau = 1

# Data generation
beta.raw = 1/sqrt(1:p)
beta.main =  C  * beta.raw / sqrt(sum(beta.raw^2))
delta.clust = 4 / sqrt(n) * rep(1, p)
gen.data = function() {
  CLUST = rbinom(n, 1, 0.5)
  W = rbinom(n, 1, eps + (1 - 2*eps) * CLUST)
  X = matrix(rnorm(n * p), n, p) + outer(CLUST, delta.clust)
  Y = X %*% beta.main + rnorm(n) + tau * W
  list(X=X, Y=Y, W=W, catt=tau)
}
data = gen.data()
X = data$X
Y = data$Y
W = data$W

# Run residual balancing
tau.rb = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "elnet", alpha = 0.9, zeta = 0.5)
tau.rb.plain = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "none", zeta = 0.5)


tau.rb
tau.rb.plain
