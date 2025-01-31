library(ebalTorchkmk)
library(glue)

# Hyperparameters
n = 500
p = 100
eps = 0.1
C = 1

tau = 10

# Data generation
beta.raw = 1/sqrt(1:p)
beta.main =  C  * beta.raw / sqrt(sum(beta.raw^2))
delta.clust = 10 / sqrt(n) * rep(1, p)
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

# Define modified EBAL
ebal_total = function(Treatment,
                      X,
                      base.weight = NULL,
                      norm.constant = NULL,
                      coefs = NULL,
                      max.iterations = 200,
                      constraint.tolerance = 1,
                      print.level = 0,
                      method = c("GaussNewton", "AutoDiff")){
  ######################################################################
  # Checks
  ######################################################################
  if (sum(Treatment != 1 & Treatment != 0) > 0) {
    stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
  }
  if (var(Treatment) == 0) {
    stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
  }
  X = as.matrix(X)
  Treatment = as.numeric(Treatment)
  if (sum(is.na(X)) > 0) stop("X contains missing data")
  if (sum(is.na(Treatment)) > 0) stop("Treatment contains missing data")
  if (length(Treatment) != nrow(X)) stop("length(Treatment) != nrow(X)")
  if (length(max.iterations) != 1) stop("length(max.iterations) != 1")
  if (length(constraint.tolerance) != 1) stop("length(constraint.tolerance) != 1")
  # set up elements
  ntreated = sum(Treatment == 1)
  ncontrols = sum(Treatment == 0)
  # Base weight is used for all samples
  if (is.null(base.weight)) base.weight = rep(1, (ncontrols + ntreated))
  if (length(base.weight) != (ncontrols + ntreated)) {
    stop("length(base.weight) !=  number of controls  sum(Treatment==0)")
  }
  
  ##################################
  ######### Call Optimizer #########
  ##################################
  
  # Support only AutoDiff for now
  if (method == 'AutoDiff'){
    target = colMeans(X)
    # Weight for the control group
    ctr_source = X[Treatment == 0, ]
    # torch fitter
    eb.out = ebalTorchkmk::ebal_torch(
      X0 = ctr_source,       # donor matrix
      X1 = target,           # target moments
      # base_weight = base.weight[Treatment == 0]/ncontrols,
      maxit = max.iterations
    )
    # return list - fewer elements than other option
    z_ctr = list(
      target.margins = target,
      co.xdata = ctr_source,
      w = eb.out$Weights.ebal,
      coefs = eb.out$coefs,
      max.iterations = max.iterations,
      base.weight = base.weight[Treatment==0]
    )
    class(z_ctr) = "ebalance"
    
    # Weight for the treatment group
    trt_source = X[Treatment == 1, ]
    # torch fitter
    eb.out = ebalTorchkmk::ebal_torch(
      X0 = trt_source,       # donor matrix
      X1 = target,           # target moments
      # base_weight = base.weight[Treatment==1]/ntreated,
      maxit = max.iterations
    )
    # return list - fewer elements than other option
    z_trt = list(
      target.margins = target,
      co.xdata = trt_source,
      w = eb.out$Weights.ebal,
      coefs = eb.out$coefs,
      max.iterations = max.iterations,
      base.weight = base.weight[Treatment==1]
    )
    class(z_trt) = "ebalance"
    
    return(list(control=z_ctr, treated=z_trt))
  }
  else {
    stop('Only AutoDiff is supported for now')
  }
}

# Run
ebal_res = ebal_total(Treatment = W, X=X, 
                          constraint.tolerance = 1, # default tolerance = 1
                          method = 'AutoDiff'
                          )
ebal_w_ctr = ebal_res$control$w
ebal_w_trt = ebal_res$treated$w

tau.ebal = sum(ebal_w_trt[W == 1]*Y[W == 1]) - sum(ebal_w_trt[W == 1] * Y[W == 1])
tau.raw = mean(Y[W == 1]) - mean(Y[W == 0])

print('Final ATE')
print(glue('raw: {tau.raw}'))
print(glue('ebal: {tau.ebal}'))