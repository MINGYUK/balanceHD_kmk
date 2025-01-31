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
NREP = 1
points = 9
k=p/10

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
} else if(beta.setup == 7){
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
  # results_list_raw[[length(results_list_raw) + 1]] = colMeans(t(results))
}

# Convert results to dataframe
results_df = do.call(rbind, results_list)
print(results_df)
write.csv(results_df, file=glue('final_df_{p}_{NREP}_{points}.csv'))

# Plot MSE vs. Zeta

library(ggplot2)
library(tidyr)

# Convert matrix to data frame and exclude CATT column
plot_data <- as.data.frame(results_df[, -1])

# Create x-axis values from your custom X list
plot_data$x <- zeta_values  # Replace X with your actual vector of coordinates

# Reshape data to long format for plotting
plot_data_long <- pivot_longer(plot_data, 
                               cols = -x, 
                               names_to = "Method", 
                               values_to = "Value")

# Create the plot with centered title and d=100 annotation
ggplot(plot_data_long, aes(x = x, y = Value, color = Method)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Tolerance / Regularization scaling", y = "Mean Squared Error", 
       title = glue("d = {p}")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "bottom"
  )
# annotate("label", x = Inf, y = Inf, label = glue("d = {p}"),  # Add annotation
#          hjust = 1.1, vjust = 1.5, color = "black", size = 4)

# Save the plot
ggsave(glue("mse_{p}_{NREP}_{points}.png"), width = 8, height = 6, dpi = 300)