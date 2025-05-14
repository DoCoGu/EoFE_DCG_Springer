# CAPM.R
# This script tests the single-factor Capital Asset Pricing Model (CAPM)
# for a set of assets by regressing excess asset returns on the market
# excess return in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Clear the workspace (optional, but good practice)
rm(list = ls())

# Load required libraries
library(readxl) # For reading Excel files

# --- Data Loading ---
# Import Stock Returns from an Excel file
Ret <- read_excel('../Data/Returns.xlsx')
# Import Fama-French Factors from an Excel file
Factors <- read_excel('../Data/FF_Factors.xlsx')

# --- Data Preparation ---
# Extract Risk-Free Rate (Rf) and Market Excess Return (Mkt)
# Select the 5th column from all rows of Factors, divide by 100, convert to matrix
Rf <- as.matrix(Factors[, 5] / 100)
# Select the 2nd column from the second row onwards of Factors, divide by 100, convert to matrix
Mkt <- as.matrix(Factors[2:276, 2] / 100) # Note: Using specific row range 2:276 as in original code

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert asset return columns (all except the first, which is date) to a matrix
# Subtract the Risk-Free Rate (Rf) from row 2 to 276 to match the length of ExRet
ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]

# Convert Mkt to a matrix (if it wasn't already)
Mkt <- as.matrix(Mkt)

# --- Testing the CAPM (One-Factor Model) ---

# Get the number of assets
mS <- ncol(ExRet)

# Initialize matrices to store regression results
# Each matrix stores the coefficient and its t-statistic for each asset
alpha_S <- matrix(NA, nrow = 2, ncol = mS)  # Alpha coefficient and its t-statistic
beta_S <- matrix(NA, nrow = 2, ncol = mS)   # Beta coefficient and its t-statistic
R2_S <- numeric(mS) # Adjusted R-squared for each asset

# Loop through each asset to perform the CAPM regression
for (i in 1:mS) {
  # Fit linear model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
  # Use lm() with ExRet[, i] as the dependent variable and Mkt as the independent variable
  mdl <- lm(ExRet[, i] ~ Mkt)

  # Extract and store the regression results (coefficients and t-statistics)
  # coef(mdl)[1] is the intercept (alpha)
  # coef(mdl)[2] is the coefficient for Mkt (beta)
  alpha_S[1, i] <- coef(mdl)[1]
  alpha_S[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1] # t-statistic for intercept
  beta_S[1, i] <- coef(mdl)[2]
  beta_S[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2] # t-statistic for Mkt

  # Store the adjusted R-squared
  R2_S[i] <- summary(mdl)$adj.r.squared
}

# --- Display and Save Results ---
# Create a data frame to display the regression results
# Each column represents an asset, and rows represent statistics (coefficient, t-stat, R2)
Results <- data.frame(
  alpha = c(alpha_S[1,]), # Alpha coefficients
  `(alpha t-stat)` = c(alpha_S[2,]), # Alpha t-statistics
  beta = c(beta_S[1,]), # Beta coefficients
  `(beta t-stat)` = c(beta_S[2,]), # Beta t-statistics
  `Adj. R2` = R2_S # Adjusted R-squared
)

# Transpose the data frame to have statistics as rows and assets as columns
Results <- t(Results)

# Define row names (statistics)
rownames(Results) <- c("alpha", "(alpha t-stat)", "beta", "(beta t-stat)", "Adj. R2")
# Define column names (asset names) from the original Excess Returns matrix
colnames(Results) <- colnames(ExRet)

# Write the results table to a CSV file
write.csv(Results, file = "CAPM_Stock.csv", row.names = TRUE)


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
