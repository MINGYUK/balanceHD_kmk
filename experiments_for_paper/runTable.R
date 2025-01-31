rm(list = ls())

library(ggplot2)
library(glue)
library(tidyr)
library(dplyr)
library(ebalTorchkmk)
library(mvtnorm)
library(ggplot2)

source("/home/manjimin/Rcodes/customCodes/balanceHD_kmk/experiments_for_paper/run.comparison.R")

beta.setup = 1
prop.setup = 4
n = 500
p = c(100, 500, 1000)
eps = 0.1
C = 1
NREP = 2
points = 9

extra.param = 1
experiment = 1

tau = 1

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
    Y = X %*% beta.main + tau * W
    Y_prime = X %*% beta.main
    list(X=X, Y=Y, W=W, catt=tau, Y_prime=Y_prime)
  }
  
} else if(prop.setup == 2) {
  
  # AR-1 model
  rho = c(0.5, 0.9)[extra.param]
  beta.main = beta.main[1 + (23 * (0:(length(beta.main) - 1))) %% length(beta.main)]
  beta.prop = c(rep(1/40, 100), rep(0, p - 100))
  gen.data = function() {
    X = 10 * rmvnorm(n, sigma = outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y))), method = "chol")
    W = rbinom(n, 1, prob = 1/(1 + exp(-beta.prop * X)))
    Y = X %*% beta.main + tau * W
    Y_prime = X %*% beta.main
    list(X=X, Y=Y, W=W, catt=tau, Y_prime=Y_prime)
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
    Y = X %*% beta.main + tau.clust[cluster] * W
    Y_prime = X %*% beta.main
    list(X=X, Y=Y, W=W, catt=mean(tau.clust[cluster[W==1]]), Y_prime=Y_prime)
  }
  
} else if(prop.setup == 4) {
  
  tau = 1
  
  # misspecified model
  gen.data = function() {
    large_n = 50000
    X = matrix(rnorm(large_n * p), large_n, p)
    # Normalized such that E[tauX | W = 1] = 1
    tauX = log(1 + exp(-2 - 2 * X[,1])) / 0.915
    ptreat = 1 - exp(-tauX)
    W = rbinom(large_n, 1, prob = ptreat)
    Y = X %*% beta.main + tauX * (2 * W - 1) / 2
    Y1 = X %*% beta.main + tauX * (2 * rep(1, nrow(X)) - 1) / 2
    Y0 = X %*% beta.main + tauX * (2 * rep(0, nrow(X)) - 1) / 2
    
    print(mean(Y1) - mean(Y0))
    print(mean(tauX))
    
    set.seed(123)  # For reproducibility
    sampled_indices <- sample(1:large_n, n)
    
    X_sampled <- X[sampled_indices, ]
    Y_sampled <- Y[sampled_indices]
    W_sampled <- W[sampled_indices]
    
    list(X=X_sampled, Y=Y_sampled, W=W_sampled, catt=(mean(Y1) -mean(Y0)))
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
    Y = X %*% beta.main + W * tau
    Y_prime = X %*% beta.main
    list(X=X, Y=Y, W=W, catt=tau)
  }
  
}

if (experiment == 2) {
  
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

# Define tolerance ranges
points = points # odd
power_of = 5
zeta_values = seq(1/(points+1), 1 - 1/(points+1), length.out = points)  # Centered at 0.5
ebal_tol_values = power_of^((1:points)-(points%/%2 + 1))  # Centered at 1
sbw_tol_values = power_of^((1:points)-(points%/%2 + 1)) / 5  # Centered at 0.2

# Store results
results_list = list()

for (idx in 1:length(zeta_values)) {
  print(glue('point {idx}'))
  ebal_tol = ebal_tol_values[[idx]]
  sbw_tol = sbw_tol_values[[idx]]
  zeta = zeta_values[[idx]]
  results = replicate(NREP, {
    XYW = gen.data()
    res = run.comparison(XYW, ebal_tol, sbw_tol, zeta)
    return(c(CATT=XYW$catt, res))
  })
  
  results_list[[length(results_list) + 1]] = colMeans(
    ((t(results) - results["CATT",])/results["CATT",])^2)
}

# Convert results to dataframe
results_df = do.call(rbind, results_list)
print(results_df)
write.csv(results_df, file='final_df_{p}_{NREP}_{points}.csv')



# # Plot MSE vs. Zeta
# 
# library(ggplot2)
# library(tidyr)
# 
# # Convert matrix to data frame and exclude CATT column
# plot_data <- as.data.frame(results_df[, -1])
# 
# # Create x-axis values from your custom X list
# plot_data$x <- zeta_values  # Replace X with your actual vector of coordinates
# 
# # Reshape data to long format for plotting
# plot_data_long <- pivot_longer(plot_data, 
#                                cols = -x, 
#                                names_to = "Method", 
#                                values_to = "Value")
# 
# 
# # Create the plot with centered title and d=100 annotation
# ggplot(plot_data_long, aes(x = x, y = Value, color = Method)) +
#   geom_point(size = 3) +
#   geom_line(linewidth = 1) +
#   scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
#   labs(x = "Tolerance / Regularization scaling", y = "Mean Squared Error",
#        title = "Comparison of Balance Metrics") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5),  # Center the title
#     legend.position = "bottom"
#   ) +
#   annotate("label", x = Inf, y = Inf, label = glue("d = {p}"),  # Add annotation
#            hjust = 1.1, vjust = 1.5, color = "black", size = 4)
# 
# # Save the plot
# ggsave(glue("mse_{p}_{NREP}_{points}.png"), width = 8, height = 6, dpi = 300)
