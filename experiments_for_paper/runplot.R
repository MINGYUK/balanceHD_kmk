rm(list = ls())

library(ggplot2)
library(glue)
library(tidyr)
library(dplyr)
library(ebalTorchkmk)
library(mvtnorm)
library(ggplot2)

source("/home/manjimin/Rcodes/customCodes/balanceHD_kmk/experiments_for_paper/run.comparison.R")

set.seed(123)

beta.setup = 7
prop.setup = 4
n = 500
p = 1000
eps = 0.1
C = 1
NREP = 5
points = 9

tau = 1
df_list = list()
zeta_values = seq(1/(points+1), 1 - 1/(points+1), length.out = points)
d_values = c(100, 500, 1000)

for(p in d_values){# Intialize beta
  
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
  } else if(beta.setup == 7){
    k = p / 2
    # Generate a k-sparse beta vector (randomly selecting k coordinates)
    beta.raw = rep(0, p)
    beta_indices = sample(1:p, k)  # Pick k random coordinates
    beta.raw[beta_indices] = 1  # Set those to 1, others remain 0
  }
  
  beta.main =  C  * beta.raw / sqrt(sum(beta.raw^2))
  
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
    
    sampled_indices <- sample(1:large_n, n)
    
    X_sampled <- X[sampled_indices, ]
    Y_sampled <- Y[sampled_indices]
    W_sampled <- W[sampled_indices]
    
    list(X=X_sampled, Y=Y_sampled, W=W_sampled, ate=(mean(Y1) -mean(Y0)))
  }
  
  # Define tolerance ranges
  points = points # odd
  power_of = 5
  zeta_values = seq(1/(points+1), 1 - 1/(points+1), length.out = points)  # Centered at 0.5
  ebal_tol_values = power_of^((1:points)-(points%/%2 + 1))  # Centered at 1
  sbw_tol_values = power_of^((1:points)-(points%/%2 + 1)) / 5  # Centered at 0.2
  
  # Store results
  results_list = list()
  results_list_raw = list()
  XYW = gen.data()
  
  for (idx in 1:length(zeta_values)) {
    print(glue('point {idx}'))
    ebal_tol = ebal_tol_values[[idx]]
    sbw_tol = sbw_tol_values[[idx]]
    zeta = zeta_values[[idx]]
    results = replicate(NREP, {
      res = run.comparison(XYW, ebal_tol, sbw_tol, zeta)
      return(c(ATE=XYW$ate, res))
    })
    
    results_list[[length(results_list) + 1]] = colMeans(
      ((t(results) - results["ATE",])/results["ATE",])^2)
  }
  
  results_df = do.call(rbind, results_list)
  write.csv(results_df, file=glue(
    '/home/manjimin/Rcodes/customCodes/ICML_outcome/final_df_{p}_{NREP}_{points}.csv'))
  df_list[[length(df_list) + 1]] = results_df
  
  }

# Plot MSE vs. Zeta

library(ggplot2)
library(tidyr)
library(patchwork)

plot_list <- list()

for (i in 1:3) {
  current_df <- df_list[[i]]
  current_d <- d_values[i]
  
  # Process each dataframe
  plot_data <- as.data.frame(current_df[, -1])  # Remove first column (CATT)
  plot_data$x <- zeta_values  # Add x-axis values
  
  plot_data_long <- pivot_longer(plot_data,
                                 cols = -x,
                                 names_to = "Method",
                                 values_to = "Value")
  
  # Create individual plot
  p <- ggplot(plot_data_long, aes(x = x, y = Value, color = Method)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05))  # Add 5% right padding
    ) +
    scale_y_continuous(limits = c(0, 2)) +
    labs(
      x = ifelse(i == 2, "Tolerance / Regularization scaling", ""),
      y = ifelse(i == 1, "Mean Squared Error", ""),
      title = glue("d = {current_d}")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  plot_list[[i]] <- p
}

# Combine plots with shared legend
combined_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save combined plot
ggsave("/home/manjimin/Rcodes/customCodes/ICML_outcome/combined_mse_plots.png", combined_plot,
       width = 14, height = 5, dpi = 300)