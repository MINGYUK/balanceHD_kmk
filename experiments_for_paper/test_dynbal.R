library(dynbal)
library(balanceHDkmk)
library(glue)

# Hyperparameters
n = 500
p = 1000
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

ebal_res = dynbal::ebalance(Treatment= W,
                            X=X,
                            lr = 0.2)

ebal_w_ctr = ebal_res$control
ebal_w_trt = ebal_res$treated

tau.dynbal = sum(
  ebal_w_trt * Y[W == 1]) - sum(
    ebal_w_ctr * Y[W == 0])
print(glue('dynbal: {tau.dynbal}'))

####################### EBAL #############################
ebal_total = function(Treatment,
                      X,
                      base.weight = NULL,
                      norm.constant = NULL,
                      coefs = NULL,
                      max.iterations = 200,
                      constraint.tolerance = 1,
                      print.level = 0,
                      method = c("GaussNewton", "AutoDiff"),
                      kappa = 1,
                      lr = 0.05){
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
    w_ctr = dynbal::ebal_torch(
      X = ctr_source,       # donor matrix
      M_hat = target,           # target moments
      maxit = max.iterations,
      kappa = kappa,
      lr=lr
    )
    
    # Weight for the treatment group
    trt_source = X[Treatment == 1, ]
    # torch fitter
    w_trt = dynbal::ebal_torch(
      X = trt_source,       # donor matrix
      M_hat = target,           # target moments
      maxit = max.iterations,
      kappa = kappa,
      lr=lr
    )
    
    return(list(control=w_ctr, treated=w_trt))
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
ebal_w_ctr = ebal_res$control
ebal_w_trt = ebal_res$treated

tau.ebal = sum(ebal_w_trt[W == 1]*Y[W == 1]) - sum(ebal_w_trt[W == 1] * Y[W == 1])

print(glue('ebal: {tau.ebal}'))
