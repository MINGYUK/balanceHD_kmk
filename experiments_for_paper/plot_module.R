rm(list = ls())

library(ebalTorchkmk)
library(mvtnorm)
library(ggplot2)

args=(commandArgs(TRUE))
beta.setup = as.numeric(args[1])
prop.setup = as.numeric(args[2])
n = as.numeric(args[3])
p = as.numeric(args[4])
eps = as.numeric(args[5])
C = as.numeric(args[6])

NREP = 50
extra.param = 1
experiment = 1

tau = 1

source("./run.comparison.R")

# Define tolerance ranges
zeta_values = seq(0.1, 0.9, by=0.1)  # Centered at 0.5
ebal_tol_values = c(0.25, 0.5, 1, 2, 4)  # Centered at 1
sbw_tol_values = c(0.005, 0.01, 0.02, 0.04, 0.08)  # Centered at 0.02

# Store results
results_list = list()

for (zeta in zeta_values) {
  for (ebal_tol in ebal_tol_values) {
    for (sbw_tol in sbw_tol_values) {
      
      mse_values = replicate(NREP, {
        XYW = gen.data()
        res = run.comparison(XYW, ebal_tol, sbw_tol, zeta)
        mse = mean(((res - XYW$catt) / XYW$catt)^2)  # Compute MSE
        return(mse)
      })
      
      avg_mse = mean(mse_values)  # Average MSE over 50 repetitions
      
      results_list = append(results_list, list(data.frame(
        Zeta = zeta,
        EBAL_Tol = ebal_tol,
        SBW_Tol = sbw_tol,
        MSE = avg_mse
      )))
    }
  }
}

# Combine results
results_df = do.call(rbind, results_list)

# Plot results
ggplot(results_df, aes(x = Zeta, y = MSE, color = factor(EBAL_Tol), shape = factor(SBW_Tol))) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = zeta_values) +
  labs(title = "MSE of EBAL, SBW, and Residual Balance",
       x = "Zeta",
       y = "Mean Squared Error (MSE)",
       color = "EBAL Tolerance",
       shape = "SBW Tolerance") +
  theme_minimal()

# Save plot
ggsave("mse_plot.png", width = 8, height = 6, dpi = 300)

save.image(file=paste0("results/simulation_results.RData"))
