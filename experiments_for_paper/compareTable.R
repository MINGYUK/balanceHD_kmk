library(balanceHDkmk)
library(ebalTorchkmk)
library(Matrix)

compare.table = function(XYW, ebal_tol, sbw_tol, resHD_zeta) {
  
  X = XYW$X
  Y = XYW$Y
  W = XYW$W
  
  # Run residual balancing
  print('running residHD')
  tau.rb = residualBalance.ate(X, Y, W, target.pop = c(0, 1), fit.method = "elnet", alpha = 0.9, zeta = resHD_zeta)
  # tau.rb.plain = residualBalance.ate(X, Y, W, target.pop = c(0, 1), fit.method = "none", zeta = resHD_zeta)
  
  ################### Custom implementation by manjimin
  # Run EBAL
  print('running ebal')
  tau.ebal = ebal.ate(Treatment = W,
                      X = X,
                      Y=Y,
                      method='AutoDiff',
                      constraint.tolerance = ebal_tol)
  print(tau.ebal)
  
  print('running sbw')
  # Run SBW
  tau.sbw = sbw.ate(W,
                    X,
                    Y,
                    tolerance = sbw_tol,
                    algorithm = FALSE,
                    std = 'group',
                    solver = 'osqp')
  
  print('running dynbal')
  # Run dynbal
  tau.dynbal = dynbal.ate(Treatment = W,
                     X = X,
                     Y=Y,
                     method='AutoDiff',
                     constraint.tolerance = 1)
  
  #####################################################
  
  results = c(
    ARB=tau.rb,
    # ResidBalancePlain=tau.rb.plain,
    EBAL = tau.ebal,
    SBW = tau.sbw,
    DYNBAL = tau.dynbal
  )
  return(results)
  
}

coverage.test = function(XYW) {
  
  X = XYW$X
  Y = XYW$Y
  W = XYW$W
  
  # Run residual balancing
  tau.rb = residualBalance.ate(X, Y, W, target.pop = 1, fit.method = "elnet", alpha = 0.9, zeta = 0.5, estimate.se = TRUE)
  return(tau.rb)
}
