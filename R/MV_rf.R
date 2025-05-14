# MV_rf.R
# This script implements Mean-Variance portfolio optimization with a
# risk-free asset, deriving and plotting the Capital Market Line (CML),
# and calculates optimal portfolio weights in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Load required libraries
library(readxl)  # For reading Excel files
library(writexl) # For writing Excel files
library(ggplot2) # For creating plots
library(dplyr)   # For data manipulation (used for colMeans and select)

# --- Data Loading ---
# Load returns data from an Excel file
Ret <- read_excel('../Data/Returns.xlsx')
# Load Fama-French Factors from an Excel file
Factors <- read_excel('../Data/FF_Factors.xlsx')

# --- Data Preparation ---
# Extract Risk-Free Rate (Rf)
# Select the 5th column from the second row onwards of Factors
# Divide by 1200 to convert from annual percentage to monthly decimal rate
Rf <- as.matrix(Factors[-1, 5] / 1200) # Use as.matrix to ensure matrix format

# Extract asset return data (assuming it starts from the 2nd column)
# Convert to a matrix
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

# --- Calculate Parameters for the Capital Market Line (CML) ---
# Calculate the mean of the risk-free rate
mean_rf <- mean(Rf)

# Calculate H, which is related to the squared Sharpe Ratio of the tangency portfolio
# H = (z - rf*1)' * inv(S) * (z - rf*1)
# z is a vector, mean(Rf) is a scalar. z - mean(Rf) is a vector.
# t(vector) %*% matrix %*% vector is correct matrix multiplication.
H <- t(z - mean(Rf)) %*% V1 %*% (z - mean(Rf))
sqrtH <- sqrt(as.numeric(H)) # Extract the scalar value from H and take the square root

# Define a range of target portfolio expected returns for plotting the CML
# Create a sequence from the mean risk-free rate up to 0.009, with 100 points
mu_p_cml <- seq(mean(Rf), 0.009, length.out = 100)

# Calculate the variance of portfolios on the Capital Market Line
# Formula: (1/H) * (mu_p - rf)^2
# mu_p_cml is a vector, mean(Rf) is scalar. (mu_p_cml - mean(Rf))^2 performs element-wise squaring.
sig2_p_cml <- 1 / as.numeric(H) * (mu_p_cml - mean(Rf))^2
# Calculate the standard deviation (risk) of portfolios on the CML
sig_p_cml <- sqrt(sig2_p_cml) # Element-wise square root

# --- Calculate Tangency Point ---
# Calculate the expected return of the Tangency Portfolio
# mu_t = rf + H / (1' * inv(S) * (z - rf*1))
# rep(1, N) creates a vector of ones. V1 %*% (z - rep(mean(Rf), N)) is a vector.
# t(rep(1, N)) %*% V1 %*% (z - rep(mean(Rf), N)) is a scalar.
denominator_t <- t(rep(1, N)) %*% V1 %*% (z - rep(mean(Rf), N))
mu_t <- mean(Rf) + as.numeric(H) / as.numeric(denominator_t)

# Calculate the standard deviation of the Tangency Portfolio
# sig_t = sqrt(H) / (1' * inv(S) * (z - rf*1))
sig_t <- sqrtH / as.numeric(denominator_t)

# --- Efficient Frontier without Risk-Free Asset (for comparison) ---
# Recalculate parameters A, B, C, D for the risky asset efficient frontier
A <- t(z) %*% V1 %*% z
B <- t(z) %*% V1 %*% rep(1, N) # rep(1, N) is a column vector of ones
C <- t(rep(1, N)) %*% V1 %*% rep(1, N) # Corrected C calculation
D <- as.numeric(A) * as.numeric(C) - as.numeric(B)^2 # Ensure scalar multiplication

# Define a range of target expected returns for the risky asset efficient frontier
mu_p_risky <- seq(0.001, 0.008, length.out = 100)

# Calculate the variance of the risky asset efficient frontier
# Formula: (1/D) * (C*mu_p^2 - 2*B*mu_p + A)
sig2_p_risky <- (1 / as.numeric(D)) * (as.numeric(C) * mu_p_risky^2 - 2 * as.numeric(B) * mu_p_risky + as.numeric(A))
# Calculate the standard deviation (risk) of the risky asset efficient frontier
sig_pp_risky <- sqrt(sig2_p_risky)

# Set standard deviations below a certain threshold to NaN for cleaner plotting
# This is often done to show only the upper (efficient) part of the frontier
sig_pp_risky[mu_p_risky < 0.005] <- NaN

# --- Plotting ---
# Create a ggplot object
p <- ggplot() +
  # Add the Capital Market Line (CML)
  geom_line(aes(x = sig_p_cml, y = mean(Rf) + sqrtH * sig_p_cml, color = 'CML'), size = 1) + # Map color for legend
  # Add the risky asset efficient frontier
  geom_path(aes(x = sig_pp_risky, y = mu_p_risky, color = 'Risky Frontier'), size = 1) + # Map color for legend
  # Add a scatter point for the Tangency Portfolio
  geom_point(aes(x = sig_t, y = mu_t, color = 'Tangency Portfolio'), size = 3) + # Map color for legend
  # Set plot title and axis labels
  labs(
    title = 'Capital Market Line',
    x = 'Portfolio Risk',
    y = 'Portfolio Expected Return',
    color = "Portfolio Type" # Legend title
  ) +
  # Apply a minimal theme
  theme_minimal() +
  # Customize legend position
  theme(legend.position = 'topleft') +
  # Set x-axis limits
  xlim(0, 0.05) +
  # Add a horizontal dashed line at the expected return of the Tangency Portfolio
  geom_hline(yintercept = mu_t, linetype = 'dashed', color = 'black', size = 1) +
  # Add the lower part of the CML (borrowing/lending at risk-free rate)
  # Need to create the lower part explicitly if not handled by the first geom_line
  # The first geom_line already plots the full line from mean(Rf) upwards.
  # Let's add the line below mean(Rf) separately if needed, or rely on the first line.
  # The formula mean(Rf) + sqrtH * sig_p_cml covers both lending (sig_p_cml > 0) and borrowing (sig_p_cml < 0, though std dev is non-negative).
  # The line mean(Rf) - sqrtH * sig_p_cml seems redundant if sqrtH is the positive root.
  # Let's assume the intent was to show the line extending below the risk-free rate.
  # The formula for the line through Rf and Tangency is E[R_p] = Rf + (E[R_T] - Rf)/sigma_T * sigma_p
  # Which is Rf + sqrt(H) * sigma_p. This is already plotted.
  # The line `aes(x = sig_p, y = mean(Rf) - sqrtH * sig_p)` seems incorrect if sqrtH is positive.
  # Let's remove the potentially incorrect lower line and rely on the first CML line.
  # geom_line(
  #   aes(x = sig_p_cml, y = mean(Rf) - sqrtH * sig_p_cml), # Removed this line
  #   color = 'blue',
  #   size = 1,
  #   linetype = 'solid'
  # ) +
  # Add a vertical dashed line at the standard deviation of the Tangency Portfolio
  geom_vline(xintercept = sig_t, linetype = 'dashed', color = 'black', size = 1) +
  # Define custom colors for the lines and points
  scale_color_manual(values = c("CML" = "blue", "Risky Frontier" = "black", "Tangency Portfolio" = "red"))


# Save the plot as a PNG file
ggsave('MV_rf.png', plot = p, device = 'png') # Specify plot object

# --- Optimal weights Calculation ---
# Define a range of target expected returns for which to calculate optimal weights on the CML
mup <- seq(0.001, 0.05, length.out = 10)

# Calculate optimal weights for portfolios on the CML (combination of risk-free and tangency)
# Formula: w_p = inv(S) * (z - rf*1) * (mu_p - rf) / H
# V1 is inv(S), (z - mean(Rf)) is (z - rf*1), H is H.
# (mup_target - mean(Rf)) / as.numeric(H) is a vector of multipliers for each target return.
# We need to multiply the vector V1 %*% as.matrix(z - mean(Rf)) by each element of this multiplier vector.
# In R, matrix %*% vector gives a vector. Matrix %*% matrix gives matrix multiplication.
# To get a matrix of weights (N x length(mup_target)), we can use element-wise multiplication
# or matrix multiplication with a diagonal matrix.
# Let's use matrix multiplication with a diagonal matrix created from the multipliers.
multipliers <- as.vector((mup_target - mean(Rf)) / as.numeric(H))
wp <- V1 %*% (z - rep(mean(Rf), N)) %*% as.numeric((mup - mean(Rf)) / H)

# Tangency portfolio weights (risky assets only)
# Formula: wT = inv(S)*(z - rf*1) / (1'*inv(S)*(z - rf*1))
# V1 %*% as.matrix(z - mean(Rf)) is V1 * (z - rf*1)
# denominator_t is 1' * inv(S) * (z - rf*1)
wT <- (V1 %*% (z - rep(mean(Rf), N))) / sum(V1 %*% (z - rep(mean(Rf), N)))

# Weight on the risk-free asset
# For each portfolio (Tangency and the ones from mup_target), the RF weight is 1 - sum of risky asset weights
# cbind(wT, wp) combines Tangency weights (1 column) and wp weights (length(mup_target) columns)
# colSums calculates the sum of weights for risky assets for each portfolio
rf_weights <- 1 - colSums(cbind(wT, wp))

# --- Create Weights table ---
# Combine all calculated weights into a single data frame
# rbind combines rows: risky asset weights (wT and wp), then risk-free weights (rf_weights)
# cbind is used within rbind to combine wT and wp side-by-side
Weights <- data.frame(rbind(cbind(wT, wp), rf_weights))

# Define column names for the weights table
# First column is "GMVP" (should be "Tangency" or similar for MV_rf), followed by target return percentages
# Let's name the first column "Tangency" as it represents the Tangency Portfolio weights
colnames(Weights) <- c("Tangency", paste0(round(mup_target * 100, 2), "%"))

# Define row names for the weights table (asset names + Rf)
rownames(Weights) <- c(colnames(Ret)[-1], "Rf") # Asset names from original Returns DataFrame, plus "Rf"

# --- Save Weights to Excel ---
# Write the weights table to an Excel file
write_xlsx(Weights, "Weights_MV_rf.xlsx")


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
