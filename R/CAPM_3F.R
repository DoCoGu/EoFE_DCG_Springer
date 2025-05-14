# CAPM_3F.R
# This script tests the Fama-French three-factor model for a set of assets
# by regressing excess asset returns on the market, SMB, and HML factors in R.
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

# --- Data Preparation ---
# Extract Fama-French factors (Mkt-RF, SMB, HML)
# Select columns 2 to 4 from the second row onwards of Factors
# Divide by 100 to convert from percentage to decimal
FF <- Factors[2:276, 2:4] / 100 # Note: Using specific row range 2:276 as in original code

# Extract Risk-Free Rate (Rf)
# Select the 5th column from all rows of Factors
# Divide by 100 to convert from percentage to decimal
Rf <- Factors[, 5] / 100

# Extract Market Excess Return (Mkt) - Note: Mkt is already in FF, but extracted separately here
# Select the 2nd column from the second row onwards of Factors
# Divide by 100 to convert from percentage to decimal
Mkt <- Factors[2:276, 2] / 100 # Note: Using specific row range 2:276 as in original code

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert asset return columns (all except the first, which is date) to a matrix
# Subtract the Risk-Free Rate (Rf) from row 2 to 276 to match the length of ExRet
ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]

# Convert Mkt and FF to matrix format
Mkt <- as.matrix(Mkt) # Convert Mkt to matrix
FF <- as.matrix(FF)   # Convert FF to matrix

# --- Testing the Fama-French 3-Factor Model ---

# Get the number of assets
mS <- ncol(ExRet)

# Initialize matrices to store regression results
# Each matrix stores the coefficient and its t-statistic for each asset
alpha <- matrix(NA, nrow = 2, ncol = mS)  # Alpha coefficient and its t-statistic
beta1 <- matrix(NA, nrow = 2, ncol = mS)  # Beta for Market factor and its t-statistic
beta2 <- matrix(NA, nrow = 2, ncol = mS)  # Beta for SMB factor and its t-statistic
beta3 <- matrix(NA, nrow = 2, ncol = mS)  # Beta for HML factor and its t-statistic
R2_S <- numeric(mS) # Adjusted R-squared for each asset

# Loop through each asset to perform the Fama-French 3-Factor regression
for (i in 1:mS) {
  # Fit linear model: ExRet_i = alpha_i + beta1_i*Mkt + beta2_i*SMB + beta3_i*HML + epsilon_i
  # Use lm() with ExRet[, i] as the dependent variable and FF as the independent variables
  mdl <- lm(ExRet[, i] ~ FF)

  # Extract and store the regression results (coefficients and t-statistics)
  # coef(mdl)[1] is the intercept (alpha)
  # coef(mdl)[2] is the coefficient for FF[,1] (Mkt)
  # coef(mdl)[3] is the coefficient for FF[,2] (SMB)
  # coef(mdl)[4] is the coefficient for FF[,3] (HML)
  alpha[1, i] <- coef(mdl)[1]
  alpha[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1] # t-statistic for intercept
  beta1[1, i] <- coef(mdl)[2]
  beta1[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2] # t-statistic for Mkt
  beta2[1, i] <- coef(mdl)[3] # Corrected: Use coef(mdl)[3] for beta2 (SMB)
  beta2[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3] # t-statistic for SMB
  beta3[1, i] <- coef(mdl)[4] # Corrected: Use coef(mdl)[4] for beta3 (HML)
  beta3[2, i] <- summary(mdl)[["coefficients"]][, "t value"][4] # t-statistic for HML

  # Store the adjusted R-squared
  R2_S[i] <- summary(mdl)$adj.r.squared
}

# --- Display and Save Results ---
# Create a data frame to display the regression results
# Each column represents an asset, and rows represent statistics (coefficient, t-stat, R2)
Stocks <- data.frame(
  alpha = c(alpha[1,]), # Alpha coefficients
  `(alpha t-stat)` = c(alpha[2,]), # Alpha t-statistics
  mkt = c(beta1[1,]), # Market Beta coefficients
  `(mkt t-stat)` = c(beta1[2,]), # Market Beta t-statistics
  smb = c(beta2[1,]), # SMB Beta coefficients # Corrected: Use beta2[1,]
  `(smb t-stat)` = c(beta2[2,]), # SMB Beta t-statistics
  hml = c(beta3[1,]), # HML Beta coefficients # Corrected: Use beta3[1,]
  `(hml t-stat)` = c(beta3[2,]), # HML Beta t-statistics
  `Adj. R2` = R2_S # Adjusted R-squared
)

# Transpose the data frame to have statistics as rows and assets as columns
Stocks <- t(Stocks)

# Define row names (statistics)
rownames(Stocks) <- c("alpha", "(alpha t-stat)", "mkt", "(mkt t-stat)",
                      "smb", "(smb t-stat)", "hml", "(hml t-stat)", "Adj. R2")
# Define column names (asset names) from the original Excess Returns matrix
colnames(Stocks) <- colnames(ExRet)

# Write the results table to a CSV file
write.csv(Stocks, file = "CAPM_3F.csv", row.names = TRUE)


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
