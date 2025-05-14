# MV.R
# This script implements the Mean-Variance portfolio optimization for a set
# of risky assets and plots the efficient frontier, also calculating
# optimal portfolio weights in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Load required libraries
library(readxl)  # For reading Excel files
library(writexl) # For writing Excel files
library(ggplot2) # For creating plots
library(dplyr)   # For data manipulation (used for colMeans)

# Clear the workspace (optional, but good practice)
rm(list = ls())

# --- Data Loading and Preparation ---
# Read returns data from an Excel file
# Use col_names = TRUE to ensure the first row is treated as column names
Ret <- read_xlsx("../Data/Returns.xlsx", col_names = TRUE)

# Extract asset return data (all columns except the first, which is assumed to be Date)
# Convert to a matrix for calculations
R <- as.matrix(Ret[, -1])
N <- ncol(R) # Number of assets

# --- Calculate Statistical Measures ---
# Calculate mean returns (expected returns) for risky assets
z <- colMeans(R) # Calculate column means to get a vector of expected returns

# Calculate standard deviations (risk) for risky assets
sig <- apply(R, 2, sd) # Calculate standard deviation for each column (asset)

# Calculate the covariance matrix of returns for risky assets
V <- cov(R)
# Calculate the inverse of the covariance matrix
V1 <- solve(V)

# --- Calculate Parameters for the Efficient Frontier ---
# Calculate the parameters A, B, C, and D from the portfolio theory
# These parameters define the shape and location of the efficient frontier
# A = z' * inv(V) * z
A <- t(z) %*% V1 %*% z
# B = 1' * inv(V) * z
B <- t(rep(1, N)) %*% V1 %*% z # rep(1, N) creates a column vector of ones
# C = 1' * inv(V) * 1
C <- t(rep(1, N)) %*% V1 %*% rep(1, N) # rep(1, N) creates a column vector of ones
# D = A*C - B^2
D <- as.numeric(A) * as.numeric(C) - as.numeric(B)^2 # Ensure scalar multiplication

# Define a range of target portfolio expected returns for plotting the efficient frontier
# Create a sequence from 0 to 0.015 with a step of 0.0001
mu_p <- seq(0, 0.015, by = 0.0001)

# Calculate the variance of the efficient frontier for the range of expected returns
# Formula: (1/D) * (C*mu_p^2 - 2*B*mu_p + A)
# Use broadcasting for element-wise operations (mu_p^2, 2*B*mu_p)
sig2_p <- (1 / as.numeric(D)) * (as.numeric(C) * mu_p^2 - 2 * as.numeric(B) * mu_p + as.numeric(A))
# Calculate the standard deviation (risk) of the efficient frontier
sig_p <- sqrt(sig2_p) # Element-wise square root

# --- Plotting the Efficient Frontier ---
# Create a data frame for the efficient frontier line
eff_frontier <- data.frame(sig_p, mu_p)
colnames(eff_frontier) <- c("Portfolio Risk", "Portfolio Expected Return")

# Create a data frame for the Global Minimum Variance Portfolio (GMVP) point
# GMVP Standard Deviation = 1/sqrt(C), GMVP Expected Return = B/C
eff_frontier_point <- data.frame(
  `Portfolio Risk` = 1 / sqrt(as.numeric(C)),
  `Portfolio Expected Return` = as.numeric(B) / as.numeric(C)
)
colnames(eff_frontier_point) <- c("Portfolio Risk", "Portfolio Expected Return")

# Create a ggplot object
p <- ggplot(eff_frontier, aes(x = `Portfolio Risk`, y = `Portfolio Expected Return`)) +
  # Add the Efficient Frontier line
  geom_path(color = "black", size = 1) +
  # Add a horizontal dashed line at the expected return of the GMVP
  geom_hline(yintercept = as.numeric(B) / as.numeric(C), linetype = "dashed", color = "black", size = 1) +
  # Add a vertical dashed line at the standard deviation of the GMVP
  geom_vline(xintercept = 1 / sqrt(as.numeric(C)), linetype = "dashed", color = "black", size = 1) +
  # Add scatter points for individual stocks (Risk vs. Return)
  # Use different colors and shapes for clarity
  geom_point(data = data.frame(sig, z), aes(x = sig, y = z),
             color = c("orange", "red", "blue", "green", "pink"), # Colors for each stock
             size = 3, shape = 21, fill = c("orange", "red", "blue", "green", "pink")) + # Shape and fill for points
  # Add a scatter point for the Global Minimum Variance Portfolio (GMVP)
  geom_point(data = eff_frontier_point, aes(x = `Portfolio Risk`, y = `Portfolio Expected Return`),
             color = "black", size = 3, shape = 21, fill = "black") + # Color, size, shape, fill for GMVP point
  # Set x-axis limits
  xlim(0, 0.1) +
  # Set plot title and axis labels
  labs(
    title = "Efficient Frontier",
    x = "Portfolio Risk",
    y = "Portfolio Expected Return"
  ) +
  # Apply a minimal theme
  theme_minimal()

# Save the plot as a PNG file
ggsave("MV.png", plot = p, width = 10, height = 6) # Specify plot object and dimensions

# --- Optimal weights Calculation ---
# Calculate parameters for portfolio weights formula on the efficient frontier
# g and h are vectors that define the linear relationship between weights and expected return
# Formula: w_p = g + h * mu_p
# Ensure matrix multiplication is used correctly
g <- as.numeric(1 / D) * (as.numeric(A) * (V1 %*% rep(1, N)) - as.numeric(B) * (V1 %*% as.matrix(z))) # Corrected z to matrix
h <- as.numeric(1 / D) * (as.numeric(C) * (V1 %*% as.matrix(z)) - as.numeric(B) * (V1 %*% rep(1, N))) # Corrected z to matrix

# Define a range of target expected returns for which to calculate optimal weights
# Create a sequence from 0.001 to 0.05 with 10 points
mup_target <- seq(0.001, 0.05, length.out = 10)

# Calculate optimal weights for the target expected returns
# Formula: w_p = g + h * mu_p
# g is a vector (N x 1), h is a vector (N x 1), mup_target is a vector (1 x length(mup_target))
# To get a matrix of weights (N x length(mup_target)), we need g replicated for each target return
# and h multiplied by each target return.
# rep(g, times = length(mup_target)) replicates the vector g.
# h %*% t(mup_target) performs matrix multiplication (N x 1) %*% (1 x length) = (N x length)
wp <- rep(g, times = length(mup_target)) + h %*% t(as.matrix(mup_target)) # Ensure mup_target is a matrix for multiplication

# --- Global minimum variance portfolio weights ---
# Calculate the weights for the Global Minimum Variance Portfolio (GMVP)
# The GMVP has an expected return of B/C
w_mvp <- g + h * as.numeric(B / C) # Use as.numeric for scalar multiplication

# --- Prepare Weights Table ---
# Combine the GMVP weights and the optimal weights for target returns into a single data frame
# cbind combines columns: GMVP weights, then wp matrix
Weights <- data.frame(cbind(w_mvp, wp))

# Define column names for the weights table
# First column is "GMVP", followed by target return percentages
colnames(Weights) <- c("GMVP", paste0(round(mup_target * 100, 2), "%"))

# Define row names for the weights table (asset names)
rownames(Weights) <- colnames(Ret)[-1] # Asset names from original Returns DataFrame

# --- Save Weights to Excel ---
# Write the weights table to an Excel file
write_xlsx(Weights, "Weights_MV.xlsx") # Specify filename


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
