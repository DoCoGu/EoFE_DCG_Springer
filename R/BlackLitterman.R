# BlackLitterman.R
# This script implements the Black-Litterman model for portfolio allocation,
# blending market equilibrium returns with investor views, and plots the
# resulting efficient frontiers and calculates optimal weights in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Clear the workspace (optional, but good practice)
rm(list = ls())

# Load required libraries
library(readxl) # For reading Excel files
library(dplyr) # For data manipulation (used for column selection)
library(ggplot2) # For creating plots
# library(PerformanceAnalytics) # Often useful for financial analysis, though not explicitly used in the provided code snippet.

# Source the helper function for efficient frontier calculations
# This assumes getuncef.R is in the same directory or accessible in the R path
source("getuncef.R")

# --- Data Loading and Preparation ---

# Import asset returns from an Excel file (assuming 5 assets)
ret <- read_excel("../Data/Returns.xlsx")

# Import Fama-French Factors from an Excel file
ff <- read_excel("../Data/FF_Factors.xlsx")

# Add Market Excess Return (`Mkt-RF`) and Risk-Free Rate (`rf`) columns to the returns data frame
# Select `Mkt-RF` and `RF` columns from the second row onwards of the FF data frame
# Divide by 100 to convert from percentage to decimal
ret$`Mkt-RF` <- ff$`Mkt-RF`[2:length(ff$`Mkt-RF`)] / 100
ret$rf <- ff$RF[2:length(ff$RF)] / 100

# Calculate excess returns for assets (asset return - risk-free rate)
# Select asset return columns (from the second column up to the column before `Mkt-RF`)
asset_columns <- names(ret)[2:(ncol(ret) - 2)]
stock_returns <- as.matrix(ret[, asset_columns]) # Matrix of stock returns

# Extract the risk-free rate column
rf <- as.matrix(ret[, ncol(ret)]) # Matrix of risk-free rates

# Calculate excess returns by subtracting the risk-free rate from each stock's return
# `rep(rf, times = ncol(stock_returns))` replicates the rf column for each stock column
exret <- stock_returns - rep(rf, times = ncol(stock_returns))

# Extract Market Excess Return column
exmkt <- as.matrix(ret[, ncol(ret) - 1])

# --- Compute the Implied Equilibrium Returns (Pi) using CAPM ---
# This step estimates CAPM beta for each asset to derive equilibrium returns

T <- nrow(exret) # Number of observations
n <- ncol(exret) # Number of assets

# Initialize matrix to store CAPM betas (intercept and slope)
# Note: The original code only stored the beta (slope) coefficient.
# A standard CAPM regression gives intercept (alpha) and slope (beta).
# Let's store both, though only beta is typically used for Pi calculation.
beta_coefs <- matrix(NA, 2, n)

# Loop through each asset to perform the CAPM regression (Excess Return ~ Market Excess Return)
for (i in 1:n) {
  # Fit linear model: exret_i = alpha_i + beta_i * exmkt + epsilon_i
  # Use lm() with the current asset's excess returns and the market excess return
  mdl <- lm(exret[, i] ~ exmkt)
  # Store the coefficients (intercept and slope)
  # coef(mdl)[1] is the intercept (alpha)
  # coef(mdl)[2] is the coefficient for exmkt (beta)
  beta_coefs[, i] <- coef(mdl)
}

# Calculate the mean excess return of assets
muExret <- colMeans(exret) # Vector of mean excess returns for each asset

# Calculate the covariance matrix of the excess returns
Sigma <- cov(exret) # Covariance matrix of asset excess returns

# Calculate the implied market equilibrium excess returns (Pi) using CAPM
# Pi = beta * mean_asset_excess_return (element-wise as in the original R code)
# Note: This calculation of Pi (beta * mean_asset_excess_return) is based on the provided original R code.
# A more standard Black-Litterman approach uses beta * mean_market_excess_return.
Pi <- beta_coefs[2, ] * muExret # Element-wise multiplication of betas by mean asset excess returns

# Uncertainty factor (Litterman and He, 1999) - often set to a small value or 1/T
tau <- 1 / T

# --- Define the Investor Views (Q and P) ---
# Q is the vector of expected returns on the views, kx1 (k is the number of views)
# Example views:
Q <- c(0.04 / 12,  # View 1: Expected return for a combination of assets (monthly)
       0.02 / 12,  # View 2: Expected return for another combination (monthly)
       0.10 / 12)  # View 3: Expected return for a single asset (monthly)

# P is the pick matrix, kxn (k views, n assets)
# Defines which assets are involved in each view and their weights
# Example P matrix corresponding to the example Q views:
# View 1: Asset 1 - Asset 3
# View 2: Asset 4 - Asset 2
# View 3: Asset 5
P <- matrix(c(1, 0, -1, 0, 0,
              0, -1, 0, 1, 0,
              0, 0, 0, 0, 1), ncol = n, byrow = TRUE) # Use n for the number of assets

# Omega: Covariance matrix of the views, kxk
# Represents the uncertainty in the investor's views. It's typically diagonal
# with variances of the view errors. A common approach is P %*% (tau * Sigma) %*% t(P)
Omega <- P %*% (tau * Sigma) %*% t(P)

# --- Blend the Equilibrium returns (Pi) with the Views (Q) ---
# Calculate the Black-Litterman covariance matrix (SigmaBL)
# Formula: [(tau*Sigma)^-1 + P'*Omega^-1*P]^-1
SigmaBL <- solve(solve(tau * Sigma) + t(P) %*% solve(Omega) %*% P)

# Calculate the Black-Litterman expected returns (muBL)
# Formula: SigmaBL %*% [(tau*Sigma)^-1*Pi + P'*inv(Omega)*Q]
# Ensure Pi and Q are column vectors for matrix multiplication
Pi_col <- as.matrix(Pi) # Convert Pi to a column matrix
Q_col <- as.matrix(Q) # Convert Q to a column matrix
muBL <- SigmaBL %*% (solve(tau * Sigma) %*% Pi_col + t(P) %*% solve(Omega) %*% Q_col)

# --- Display and Save Results (Expected Returns Comparison) ---
# Create a data frame to compare mean historical excess returns, implied equilibrium returns (Pi), and Black-Litterman returns (muBL)
# Using rbind and data.frame as in the original script
muTable <- data.frame(rbind(muExret, Pi, t(muBL)))
# Set row names as in the original script
rownames(muTable) <- c("exret", "Pi", "muBL")

# Write the table to a CSV file as in the original script
write.csv(muTable, "muTable.csv", row.names = TRUE) # Write row names

# --- Implement the two efficient frontiers (Mean-Variance and Black-Litterman) ---
# Define a range of target expected returns for plotting the efficient frontiers
muR <- seq(0, 0.015, by = 0.00001)

# Calculate the Mean-Variance efficient frontier using the mean historical excess returns (muExret) and covariance (Sigma)
# Assumes getuncef() function is defined and available, returning a list with OptSigma and w
# Pass muExret as a matrix as in the original script's call to getuncef
result_MV <- getuncef(as.matrix(muExret), Sigma, muR)

# Calculate the Black-Litterman efficient frontier using the blended returns (muBL) and adjusted covariance (Sigma + SigmaBL)
# Pass muBL as a matrix as in the original script's call to getuncef
result_BL <- getuncef(muBL, Sigma + SigmaBL, muR)

# --- Create Efficient Frontier Table ---
# Create a data frame to store the efficient frontier data (returns and standard deviations)
# Using data.frame as in the original script
ef_table <- data.frame(
  r = muR,        # Target returns
  sigma = sqrt(result_MV$OptSigma), # Standard deviation for M-V frontier
  sigmaBL = sqrt(result_BL$OptSigma) # Standard deviation for B-L frontier
)

# --- Plot the Efficient Frontiers ---
# Load ggplot2 library (already loaded at the beginning, but good to ensure)
library(ggplot2)

# Create the ggplot object for the efficient frontiers plot as in the original script
EF_plot <- ggplot(ef_table, aes(x = sigma, y = r)) +
  # Add the Mean-Variance efficient frontier line with specified color and size
  geom_path(color = 'blue', size = 1) +
  # Add the Black-Litterman efficient frontier line with specified color and size
  geom_path(aes(x = sigmaBL), color = 'red', size = 1) +
  # Set plot title
  ggtitle('Efficient Frontier') +
  # Set x-axis limits
  xlim(0, 0.1) +
  # Set axis labels (using labs for both axis titles and color legend title)
  labs(x = 'Portfolio Risk', y = 'Portfolio Expected Return', color = "Frontier Type") + # Added color legend title
  # Apply a minimal theme
  theme_minimal() +
  # Customize text sizes for title and axis labels
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  # Define custom colors for the lines (mapped by aes(color = ...))
  scale_color_manual(values = c("M-V" = "blue", "B-L" = "red")) +
  # Set legend position
  theme(legend.position = "northwest") # Using "northwest" as in the original script

# Save the efficient frontier plot as an EPS file using ggsave as in the original script
ggsave("BL.eps", plot = EF_plot, device = "eps")

# --- Optimal Weights Calculation ---
# Extract the weight matrices from the getuncef results
w <- result_MV$w # M-V weights (matrix)
wBL <- result_BL$w # B-L weights (matrix)

# Define a specific target return for which to calculate optimal weights
target_return <- 0.0060

# Find the column index in the weight matrices corresponding to the target return
# `which()` returns the index where the condition is TRUE
idx <- which(ef_table$r == target_return)

# Check if the target return was found
if (length(idx) == 0) {
  print(paste("Target return", target_return, "not found in the efficient frontier range."))
} else {
  # Create a data frame to store the optimal weights for the target return as in the original script
  weights <- data.frame(matrix(ncol = n, nrow = 2))
  # Set column names to asset names as in the original script
  colnames(weights) <- colnames(ret)[2:(ncol(ret) - 2)] # Use asset_columns or derive from ret as in original
  # Set row names to identify the method (Mean-Variance or Black-Litterman)
  rownames(weights) <- c("M-V", "B-L")
  
  # Extract the weights for the target return from the weight matrices
  # w and wBL are expected to have shape (number of assets, number of target returns)
  weights[1, ] <- w[, idx] # M-V weights at target return
  weights[2, ] <- wBL[, idx] # B-L weights at target return
  
  # Write the weights table to a CSV file (using .xlsx extension as in original)
  write.csv(weights, "Weights_BL.xlsx", row.names = TRUE) # Write row names
}


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
