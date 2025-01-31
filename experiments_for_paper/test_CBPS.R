library(CBPS)

# Hyperparameters
n = 50
p = 60
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
X.mis = data$X
outcome = data$Y
treat = data$W

data_sbw = as.data.frame(cbind(W, X, Y))
names(data_sbw) = c('W', paste0("X", 1:ncol(X)), 'Y')

# CBPS
# Get PS
cbps_model = hdCBPS(treat ~ X.mis, y=outcome)
summary(cbps_model)

# Weighted regression using CBPS weights
outcome_model <- glm(Y ~ W + X, weights = cbps_model$weights)

# Summarize the results
summary(outcome_model)