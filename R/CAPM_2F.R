# CAPM_2F.R
# This script tests two-factor asset pricing models by adding either VIX
# changes or EPU changes as a second factor to the CAPM in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Clear the workspace
rm(list = ls())

# Load required libraries
library(readxl) # For reading Excel files

# --- Data Loading ---
# Import Stock Returns from an Excel file
Ret <- read_excel('../Data/Returns.xlsx')
# Import Fama-French Factors from an Excel file
Factors <- read_excel('../Data/FF_Factors.xlsx')

# Extract Risk-Free Rate (Rf) and Market Excess Return (Mkt)
# Select the 5th column from all rows of Factors, divide by 100, convert to matrix
Rf <- as.matrix(Factors[, 5] / 100)
# Select the 2nd column from the second row onwards of Factors, divide by 100, convert to matrix
Mkt <- as.matrix(Factors[2:276, 2] / 100) # Note: Using specific row range 2:276 as in original code

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert asset return columns (all except the first, which is date) to a matrix
# Subtract the Risk-Free Rate (Rf) from row 2 to 276 to match the length of ExRet
ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]

# Import Uncertainty data (including VIX and EPU) from an Excel file
Uncertainty <- read_excel('../Data/Uncertainty.xlsx')

# --- Data Preparation ---
# Calculate percentage change in VIX and EPU
# VIX change: (Current VIX / Previous VIX) - 1
# Use as.matrix to ensure matrix format, [-1,] excludes the last row, [-nrow(),] excludes the first row
VIX <- as.matrix((Uncertainty[-1, 2]) / (Uncertainty[-nrow(Uncertainty), 2]) - 1)
# EPU change: (Current EPU / Previous EPU) - 1
EPU <- as.matrix((Uncertainty[-1, 3]) / (Uncertainty[-nrow(Uncertainty), 3]) - 1)

# --- Testing the CAPM with VIX (Two-Factor Model) ---

# Get the number of assets
mS <- ncol(ExRet)

# Initialize matrices to store regression results for CAPM + VIX model
# Each matrix stores the coefficient and its t-statistic for each asset
alpha_VIX <- matrix(NA, nrow = 2, ncol = mS)  # Alpha coefficient and its t-statistic
beta1_VIX <- matrix(NA, nrow = 2, ncol = mS)  # Beta for Market factor and its t-statistic
beta2_VIX <- matrix(NA, nrow = 2, ncol = mS)  # Beta for VIX factor and its t-statistic
R2_VIX <- numeric(mS) # Adjusted R-squared for each asset

# Loop through each asset to perform the regression with Market and VIX
for (i in 1:mS) {
  # Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * VIX + epsilon_i
  # Use lm() with ExRet[, i] as the dependent variable and Mkt and VIX as independent variables
  mdl <- lm(ExRet[, i] ~ Mkt + VIX)

  # Extract and store the regression results (coefficients and t-statistics)
  # coef(mdl)[1] is the intercept (alpha)
  # coef(mdl)[2] is the coefficient for Mkt (beta1)
  # coef(mdl)[3] is the coefficient for VIX (beta2)
  alpha_VIX[1, i] <- coef(mdl)[1]
  alpha_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1] # t-statistic for intercept
  beta1_VIX[1, i] <- coef(mdl)[2]
  beta1_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2] # t-statistic for Mkt
  beta2_VIX[1, i] <- coef(mdl)[3]
  beta2_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3] # t-statistic for VIX

  # Store the adjusted R-squared
  R2_VIX[i] <- summary(mdl)$adj.r.squared
}

# --- Testing the CAPM with EPU (Two-Factor Model) ---

# Initialize matrices to store regression results for CAPM + EPU model
alpha_EPU <- matrix(NA, nrow = 2, ncol = mS)  # Alpha coefficient and its t-statistic
beta1_EPU <- matrix(NA, nrow = 2, ncol = mS)  # Beta for Market factor and its t-statistic
beta2_EPU <- matrix(NA, nrow = 2, ncol = mS)  # Beta for EPU factor and its t-statistic
R2_EPU <- numeric(mS) # Adjusted R-squared for each asset

# Loop through each asset to perform the regression with Market and EPU
for (i in 1:mS) {
  # Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * EPU + epsilon_i
  # Use lm() with ExRet[, i] as the dependent variable and Mkt and EPU as independent variables
  mdl <- lm(ExRet[, i] ~ Mkt + EPU)

  # Extract and store the regression results (coefficients and t-statistics)
  # coef(mdl)[1] is the intercept (alpha)
  # coef(mdl)[2] is the coefficient for Mkt (beta1)
  # coef(mdl)[3] is the coefficient for EPU (beta2)
  alpha_EPU[1, i] <- coef(mdl)[1]
  alpha_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1] # t-statistic for intercept
  beta1_EPU[1, i] <- coef(mdl)[2]
  beta1_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2] # t-statistic for Mkt
  beta2_EPU[1, i] <- coef(mdl)[3]
  beta2_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3] # t-statistic for EPU

  # Store the adjusted R-squared
  R2_EPU[i] <- summary(mdl)$adj.r.squared
}

# --- Display and Save Results (CAPM + VIX) ---
# Create a data frame to display the regression results for the VIX model
# Each column represents an asset, and rows represent statistics (coefficient, t-stat, R2)
Stock_VIX <- data.frame(
  alpha = c(alpha_VIX[1,]), # Alpha coefficients
  `(alpha t-stat)` = c(alpha_VIX[2,]), # Alpha t-statistics
  Mkt = c(beta1_VIX[1,]), # Market Beta coefficients
  `(Mkt t-stat)` = c(beta1_VIX[2,]), # Market Beta t-statistics
  VIX = c(beta2_VIX[1,]), # VIX Beta coefficients
  `(VIX t-stat)` = c(beta2_VIX[2,]), # VIX Beta t-statistics
  `Adj. R2` = R2_VIX # Adjusted R-squared
)

# Transpose the data frame to have statistics as rows and assets as columns
Stock_VIX <- t(Stock_VIX)

# Define row names (statistics)
rownames(Stock_VIX) <- c("alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "VIX", "(VIX t-stat)", "Adj. R2")
# Define column names (asset names) from the original Excess Returns matrix
colnames(Stock_VIX) <- colnames(ExRet)

# --- Display and Save Results (CAPM + EPU) ---
# Create a data frame to display the regression results for the EPU model
# Each column represents an asset, and rows represent statistics (coefficient, t-stat, R2)
Stock_EPU <- data.frame(
  alpha = c(alpha_EPU[1,]), # Alpha coefficients
  `(alpha t-stat)` = c(alpha_EPU[2,]), # Alpha t-statistics
  Mkt = c(beta1_EPU[1,]), # Market Beta coefficients
  `(Mkt t-stat)` = c(beta1_EPU[2,]), # Market Beta t-statistics
  EPU = c(beta2_EPU[1,]), # EPU Beta coefficients
  `(EPU t-stat)` = c(beta2_EPU[2,]), # EPU Beta t-statistics
  `Adj. R2` = R2_EPU # Adjusted R-squared
)

# Transpose the data frame to have statistics as rows and assets as columns
Stock_EPU <- t(Stock_EPU)

# Define row names (statistics)
rownames(Stock_EPU) <- c("alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "EPU", "(EPU t-stat)", "Adj. R2")
# Define column names (asset names) from the original Excess Returns matrix
colnames(Stock_EPU) <- colnames(ExRet)

# --- Write Results to CSV ---
# Write the VIX results table to a CSV file
write.csv(Stock_VIX, file = "VIX_Stock.csv", row.names = TRUE)
# Write the EPU results table to a CSV file
write.csv(Stock_EPU, file = "EPU_Stock.csv", row.names = TRUE)


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
