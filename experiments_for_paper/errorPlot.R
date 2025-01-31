rm(list = ls())

library(ggplot2)
library(tidyr)
library(patchwork)
library(glue)
library(dplyr)
library(ebalTorchkmk)
library(mvtnorm)

source("/home/manjimin/Rcodes/customCodes/balanceHD_kmk/experiments_for_paper/run.comparison.R")

set.seed(123)

beta.setup = 1
prop.setup = 4
n = 500
eps = 0.1
C = 1
NREP = 5

tau = 1
df_list = list()
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
    k = p / 10
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
  
  XYW = gen.data()
  
  # Fixed parameter settings
  ebal_tol_fixed <- 1
  sbw_tol_fixed <- 0.2
  zeta_fixed <- 0.5
  
  # Collect all repetitions
  all_errors <- replicate(NREP, {
    res <- run.comparison(XYW, ebal_tol_fixed, sbw_tol_fixed, zeta_fixed)
    relative_error <- ((res - XYW$ate)/XYW$ate)^2
    return(relative_error)
  }) %>% 
    t() %>% 
    as.data.frame() %>% 
    pivot_longer(cols = everything(),
                 names_to = "Method",
                 values_to = "MSE") %>% 
    mutate(Method = factor(Method))
  
  df_list[[length(df_list) + 1]] = all_errors
  write.table(all_errors, glue('/home/manjimin/Rcodes/customCodes/ICML_outcome/errorplot_{p}_{NREP}.csv'))
}

# Plot MSE vs. Zeta

box_plot_list <- list()
hist_plot_list = list()
for (i in 1:length(df_list)){
  current_df <- df_list[[i]]
  current_d <- d_values[[i]]
  
  # Create individual plot
    ######################
    # gather box plot
    box_plot <- ggplot(current_df, aes(x = Method, y = MSE, fill = Method)) +
    geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
    # geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
    coord_cartesian(ylim = c(0, 1.5)) +
    labs(title = glue("d = {current_d}"),
         # subtitle = glue("EBAL tolerance = {ebal_tol_fixed}, SBW tolerance = {sbw_tol_fixed}, ζ = {zeta_fixed}"),
         y = "Relative Error Squared") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ##################
  # gather hist plot
  hist_plot <-  ggplot(current_df, aes(x = MSE, fill = Method)) +
    geom_histogram(binwidth = 0.05, alpha = 0.8, position = "identity") +
    facet_wrap(~Method, nrow = 1) +
    scale_x_continuous(limits = c(0.05, 0.8))+
    coord_cartesian(ylim = c(0, NREP)) +
    labs(x = "Mean Squared Error", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  box_plot_list[[i]] <- box_plot
  hist_plot_list[[i]] <- hist_plot
}

# Combine plots with shared legend
box_combined_plot <- wrap_plots(box_plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save combined plot
ggsave(glue("/home/manjimin/Rcodes/customCodes/ICML_outcome/box_error_plot_{NREP}.png"), box_combined_plot,
       width = 14, height = 5, dpi = 300)

hist_combined_plot <- wrap_plots(hist_plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
        
ggsave(glue("/home/manjimin/Rcodes/customCodes/ICML_outcome/hist_error_plot_{NREP}.png"), hist_combined_plot,
               width = 14, height = 5, dpi = 300)