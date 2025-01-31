library(sbw)

sbw.ate = function(
    Treatment_i,
    X,
    Y,
    tolerance = 0.02,
    algorithm = FALSE,
    std = 'group',
    solver = 'osqp'
){
  data_sbw = as.data.frame(cbind(Treatment_i, X, Y))
  names(data_sbw) = c('t_ind', paste0("X", 1:ncol(X)), 'Y')
  t_ind = 't_ind'
  bal = list()
  bal$bal_cov = c(paste0("X", 1:ncol(X)))
  bal$bal_tol = tolerance # Tolerance default = 0.02
  bal$bal_std = std
  bal$bal_alg = algorithm
  
  sbw_res = sbw::sbw(dat=data_sbw, ind=t_ind, out='Y', bal=bal, 
                     sol=list(sol_nam=solver),
                     par=list(par_est='ate', par_tar=NULL))
  
  sbw.cau.est_only = function(object, out = NULL, digits, ...) {
    if (!is(object, "sbwcau")) {
      warning("Object not of class \"sbwcau\"")
      return(invisible(NULL))
    }
    ind = object$ind
    if (is.null(out)) {
      out = object$out
    }
    dat = object$dat_weights
    if (sum(1 - is.na(match(out, colnames(dat)))) == 0) {
      stop("Please specify a correct string for out.")
    }
    fac_ind = sapply(dat, is.factor)
    dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
    
    # Get ATE
    tre_ind = dat[, ind]
    weights0 = dat$sbw_weights*(1 - tre_ind)
    weights1 = dat$sbw_weights*tre_ind
    n = length(weights0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$sbw_weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    
    estimates = crossprod(weights1 - weights0, Y)
    estimates = as.vector(estimates)
    
    return(estimates)
  }
  
  sbw_est = sbw.cau.est_only(sbw_res)
  return(sbw_est)
}
