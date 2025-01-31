library(ebalTorchkmk)

# Run
ebal.ate = function(
    Treatment,
    X,
    Y,
    method,
    constraint.tolerance = 1,
    maxit = 500,
    lr = 0.2
){
  ebal_res = ebalTorchkmk::ebalance(Treatment=Treatment,
                        X=X,
                        constraint.tolerance = constraint.tolerance,
                        method=method,
                        max.iterations = maxit,
                        lr=lr)
  ebal_w_ctr = ebal_res$control$Weights.ebal
  ebal_w_trt = ebal_res$treated$Weights.ebal
  
  tau.ebal = sum(ebal_w_trt*Y[Treatment == 1]) - sum(ebal_w_ctr * Y[Treatment == 0])
  return(tau.ebal)
}