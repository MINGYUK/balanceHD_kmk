library(sbw)
# Hyperparameters
n = 500
p = 5000
eps = 0.1
C = 1

tau = 10

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

# SBW
data_sbw = as.data.frame(cbind(W, X, Y))
names(data_sbw) = c('t_ind', paste0("X", 1:ncol(X)), 'Y')
t_ind = 't_ind'
bal = list()
bal$bal_cov = c(paste0("X", 1:ncol(X)))
bal$bal_tol = 0.02 # Tolerance default = 0.02
bal$bal_std = 'group'
bal$bal_alg = FALSE

sbw_res = sbw::sbw(dat=data_sbw, ind=t_ind, out='Y', bal=bal, 
                   sol=list(sol_nam='osqp'),
                   par=list(par_est='ate', par_tar=NULL))

# Estimate output from sbwcau
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
  
  if (object$par$par_est == "att") {
    tre_ind = dat[, ind]
    weights0 = dat$sbw_weights*(1 - tre_ind)
    weights1 = dat$sbw_weights*tre_ind
    n = sum(tre_ind == 1)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$sbw_weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    
    # remove collinear columns
    qr_dat = qr(dat, tol = 1e-9, LAPACK = FALSE)
    rank_dat = qr_dat$rank
    keep_dat = qr_dat$pivot[seq_len(rank_dat)]
    dat = dat[, keep_dat, drop=FALSE]
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans(as.matrix((as.matrix(n*weights1*Y - rep(1, length(weights1))%*%(weights1%*%Y) 
                                            - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))[tre_ind == 1,])^2) + 
      colMeans(as.matrix((as.matrix(n*weights0*Y - rep(1, length(weights0))%*%(weights0%*%Y)
                                    - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))[tre_ind == 1,])^2)
    sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "atc") {
    tre_ind = dat[, ind]
    weights0 = dat$sbw_weights*(1 - tre_ind)
    weights1 = dat$sbw_weights*tre_ind
    n = sum(tre_ind == 0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$sbw_weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    
    # remove collinear columns
    qr_dat = qr(dat, tol = 1e-9, LAPACK = FALSE)
    rank_dat = qr_dat$rank
    keep_dat = qr_dat$pivot[seq_len(rank_dat)]
    dat = dat[, keep_dat, drop=FALSE]
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans(as.matrix((as.matrix(n*weights1*Y - rep(1, length(weights1))%*%(weights1%*%Y) 
                                            - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))[tre_ind == 0,])^2) + 
      colMeans(as.matrix((as.matrix(n*weights0*Y - rep(1, length(weights0))%*%(weights0%*%Y) 
                                    - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))[tre_ind == 0,])^2)
    sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "ate") {
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
    
    # remove collinear columns
    #qr_dat = qr(dat, tol = 1e-9, LAPACK = FALSE)
    #rank_dat = qr_dat$rank
    #keep_dat = qr_dat$pivot[seq_len(rank_dat)]
    #dat = dat[, keep_dat, drop=FALSE]
    #tmp = cor(dat)
    #tmp[upper.tri(tmp)] = 0
    #diag(tmp) = 0
    #dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    #
    #var_cau = colMeans((as.matrix(n*weights1*Y - rep(1, length(weights1))%*%(weights1%*%Y) 
    #                              - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))^2) + 
    #  colMeans((as.matrix(n*weights0*Y - rep(1, length(weights0))%*%(weights0%*%Y) 
    #                      - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))^2)
    #sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "cate") {
    if (is(object$par$par_tar, "character")) {
      dat = subset(dat, eval(parse(text = object$par$par_tar)))
      dat[fac_ind] = NULL
    } else if (is(object$par$par_tar, "numeric")) {
      if (sum(fac_ind) >= 1) {
        dat = dat[apply(dat[fac_ind] == object$par$par_tar[match(names(fac_ind), names(object$par$par_tar))][fac_ind], 1, prod, na.rm = TRUE) %in% 1,]
      }
      dat[fac_ind] = NULL
    }
    tre_ind = dat[, ind]
    weights0 = dat$sbw_weights*(1 - tre_ind)
    weights1 = dat$sbw_weights*tre_ind
    n = length(weights0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$sbw_weights = NULL
    dat = dat[, object$bal$bal_cov[!(object$bal$bal_cov %in% names(which(fac_ind == TRUE)))]]
    dat = as.matrix(dat)
    
    # remove collinear columns
    qr_dat = qr(dat, tol = 1e-9, LAPACK = FALSE)
    rank_dat = qr_dat$rank
    keep_dat = qr_dat$pivot[seq_len(rank_dat)]
    dat = dat[, keep_dat, drop=FALSE]
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans((as.matrix(n*weights1*Y - rep(1, length(weights1))%*%(weights1%*%Y) 
                                  - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))^2) + 
      colMeans((as.matrix(n*weights0*Y - rep(1, length(weights0))%*%(weights0%*%Y) 
                          - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))^2)
    sd_cau = sqrt(var_cau/n)
  }
  
  estimates = crossprod(weights1 - weights0, Y)
  estimates = as.vector(estimates)
  
  # cau_table = cbind(estimates, sd_cau, estimates/sd_cau, 
  #                   2*pt(q = abs(estimates/sd_cau), df = n - 1, lower.tail = FALSE),
  #                   estimates + sd_cau*qt(0.025, n-1),
  #                   estimates + sd_cau*qt(0.975, n-1))
  # rownames(cau_table) = out
  # colnames(cau_table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "2.5 %", "97.5 %")
  # cat("\n")
  # cat(paste(toupper(object$par$par_est), ": ", sep = ""), "\n")
  # print(cau_table, digits = digits)
  # cat("\n")
  # invisible(list(cau_table = cau_table))
  return(estimates)
}


sbw_est = sbw.cau.est_only(sbw_res)
sbw_est
# tau.sbw = sbw_est$cau_table[[1]]

# Custom implementation for tau
# treated = data_sbw$t_ind == 1
# control = data_sbw$t_ind == 0
# col_cov = c(paste0("X", 1:ncol(X)))
# x_t = data_sbw[,col_cov][treated,] * sbw_res$dat_weights$sbw_weights[treated]
# x_c = data_sbw[,col_cov][control,] * sbw_res$dat_weights$sbw_weights[control]
# y_t = data_sbw[,'Y'][treated]
# y_c = data_sbw[,'Y'][control]


# weighted_model <- lm(data_sbw[,'Y'] ~ data_sbw$t_ind==1, data = data_sbw[,col_cov], weights = sbw_res$dat_weights$sbw_weights)
# weighted_model