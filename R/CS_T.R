# CS_T.R
# This script performs Fama-MacBeth cross-sectional regression analysis,
# including temperature changes and consumption growth as factors,
# and calculates Newey-West and Shanken standard errors in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Load required libraries
library(readxl)    # For reading Excel files
library(sandwich)  # For robust standard errors (like Newey-West)
library(lmtest)    # For testing linear regression models and using coeftest

# Set working directory to the script's source directory
# This helps in locating data files relative to the script
setwd(getSrcDirectory(function(){})[1])

# --- Data Loading ---
# Read Portfolio Returns from an Excel file (assuming 'PortfoliosLong.xlsx')
Ret <- read_xlsx('../Data/PortfoliosLong.xlsx')
# Read Fama-French Factors from an Excel file (assuming 'FF_FactorsLong.xlsx')
FFactors <- read_xlsx('../Data/FF_FactorsLong.xlsx')
# Read Temperature data from an Excel file
Temperature <- read_xlsx('../Data/Temperature.xlsx')
# Read Consumption data from an Excel file (assuming 'ConsumptionLong.xlsx')
Consumption <- read_xlsx('../Data/ConsumptionLong.xlsx')

# --- Data Preparation ---
# Compute Temperature Change and Consumption Growth
# Temperature change: difference in log temperature
# Convert the temperature column (assuming it's the 2nd column) to a matrix and then calculate the difference in log
Temp <- as.matrix(diff(log(as.matrix(Temperature[, 2]))))
# Consumption change: percentage change (Current Consumption / Previous Consumption) - 1
# Use as.matrix to ensure matrix format, [-1,] excludes the last row, [-nrow(),] excludes the first row
C <- as.matrix(Consumption[-1, 2] / Consumption[-nrow(Consumption), 2] - 1)

# Extract Market Excess Return (Mkt) and Risk-Free Rate (Rf)
# Select the 2nd column (Mkt-RF) from the second row onwards of FFactors, convert to matrix, divide by 100
Mkt <- as.matrix(FFactors[-1, 2] / 100)
# Select the 5th column (RF) from the second row onwards of FFactors, convert to matrix, divide by 100
Rf <- as.matrix(FFactors[-1, 5] / 100)

# Calculate Excess Returns for portfolios (Portfolio Return - Risk-Free Rate)
# Select portfolio returns (all columns except the first, which is date) from the second row onwards of Ret
# Divide by 100 to convert from percentage to decimal, then subtract the Risk-Free Rate
ExRet <- as.matrix(Ret[-1, -1] / 100 - Rf)

# Combine the factors (Market, Consumption Change, Temperature Change) into a single matrix
Factors <- cbind(Mkt, C, Temp)

# --- First Stage Regression (Time-Series Regressions) ---
# Get the size of the Excess Returns and Factors matrices
n1 <- nrow(ExRet) # n1 = number of observations
n2 <- ncol(ExRet) # n2 = number of portfolios
nF <- ncol(Factors) # nF = number of factors

# Initialize matrices to store results from the first stage
# CoefAll stores the factor betas for each portfolio (excluding intercept)
CoefAll <- matrix(NA, nrow = nF, ncol = n2)
# Res stores the residuals from each regression
Res <- matrix(NA, nrow = n1, ncol = n2)

# Loop through each portfolio to perform time-series regression
for (i in 1:n2) {
  # Fit linear model: ExRet_i = alpha_i + beta_Mkt*Mkt + beta_dC*dC + beta_T*T + epsilon_i
  # Use lm() with ExRet[, i] as the dependent variable and the combined Factors matrix as independent variables
  lm_result <- lm(ExRet[, i] ~ Mkt + C + Temp)

  # Store the factor betas (excluding the intercept, which is the first coefficient)
  CoefAll[, i] <- coef(lm_result)[-1]
  # Store the residuals from the regression
  Res[, i] <- residuals(lm_result)
}

# Calculate the variance-covariance matrix of the residuals
VarCovErr <- cov(Res)

# --- Second Stage Regression (Cross-Sectional Regressions) ---
# Calculate the average excess return for each portfolio
MeanRet <- colMeans(ExRet)

# Transpose the factor betas matrix from the first stage
Betas <- t(CoefAll) # Shape: (number of portfolios, number of factors)

# Fit linear model for MeanRet against Betas
# Model: MeanRet_i = lambda_Mkt*beta_Mkt_i + lambda_dC*beta_dC_i + lambda_T*beta_T_i + eta_i
# Use lm() with `MeanRet ~ Betas - 1` to suppress the intercept (assuming lambda_0 = 0)
mdl <- lm(MeanRet ~ Betas - 1)

# Extract factor risk premia (Lambdas) and their standard errors
SE <- summary(mdl)$coefficients[, 2] # Standard Errors of the Lambdas
Lambda <- summary(mdl)$coefficients[, 1] # Factor Risk Premia (Lambdas)
Tstat <- Lambda / SE # T-statistics (standard OLS T-statistics)

# --- Shanken Correction ---
# Calculate the covariance matrix of the factors
Sigma_f <- cov(Factors)
# Calculate the inverse of (Betas' * Betas)
BtB_inv <- solve(t(Betas) %*% Betas)
# Calculate the correction term for Shanken standard errors
correction <- as.numeric(1 + t(Lambda) %*% solve(Sigma_f) %*% Lambda) # Ensure scalar result
# Calculate the Shanken-corrected variance-covariance matrix of Lambdas
# Formula: (BtB_inv * Betas' * VarCovErr * Betas * BtB_inv * correction + Sigma_f) / n1
VarLam <- (BtB_inv %*% t(Betas) %*% VarCovErr %*% Betas %*% BtB_inv * correction + Sigma_f) / n1
# Extract Shanken-corrected standard errors (square root of the diagonal elements of VarLam)
SE_Shanken <- sqrt(diag(VarLam))
# Calculate Shanken-corrected T-statistics
Tstat_Shanken <- Lambda / SE_Shanken

# --- Time-series of cross-section regressions (for Newey-West standard errors) ---
# Initialize matrix to store Lambdas from each cross-sectional regression
LambdaFull <- matrix(NA, n1, nF) # Shape: (number of observations, number of factors)

# Loop through each time period (observation) to run a cross-sectional regression
for (j in 1:n1) {
  # Use returns from a single time period for the cross-sectional regression
  MeanRet_j <- ExRet[j,]
  # Fit cross-sectional model for the current time period (using betas from the first stage)
  # Model: ExRet_j_i = lambda_Mkt_j*beta_Mkt_i + lambda_dC_j*beta_dC_i + lambda_T_j*beta_T_i + eta_j_i
  # Use lm() with `MeanRet_j ~ Betas - 1` to suppress the intercept
  mdl_j <- lm(MeanRet_j ~ Betas - 1)
  # Store the estimated Lambdas for this time period
  LambdaFull[j,] <- coef(mdl_j)
}

# Calculate the mean of the estimated Lambdas across all time periods
LambdaMean <- colMeans(LambdaFull)

# Calculate HAC-corrected (Newey-West) standard errors
# This requires the time series of estimated Lambdas (LambdaFull)
# Create a dummy independent variable (intercept only) for the Newey-West function
X_dummy <- rep(1, n1)
hac_cov <- numeric(nF) # Initialize array for HAC variances

# Loop through each factor to calculate its HAC standard error using the NeweyWest function from sandwich package
for (k in 1:nF) {
  y <- LambdaFull[, k] # Time series of estimated Lambdas for factor k
  # Fit a dummy regression (y ~ 1) to use the NeweyWest function
  mdl_dummy <- lm(y ~ X_dummy - 1) # -1 to remove the default intercept if X_dummy is just 1s
  # Compute HAC robust variance using NeweyWest
  # lag = 1 specifies the number of lags, sandwich = T for robust variance, prewhite = F, adjust = T for small sample adjustment
  hac_cov[k] <- NeweyWest(mdl_dummy, lag = 1, sandwich = TRUE, prewhite = FALSE, adjust = TRUE)[1, 1] # Extract the variance (top-left element)
}

# Calculate Newey-West Standard Errors (square root of HAC variances)
SE_NW <- sqrt(hac_cov)
# Calculate Newey-West T-statistics
Tstat_NW <- LambdaMean / SE_NW

# --- Results Table Formatting ---
# Get the names of the portfolios from the original Returns DataFrame
NamePort <- colnames(Ret)[-1]

# Create a data frame for the First Stage results (Factor Betas)
FirstStageReg <- data.frame(Mkt = Betas[, 1], dC = Betas[, 2], T = Betas[, 3])
# Set row names to portfolio names
rownames(FirstStageReg) <- NamePort

# Create a data frame for the Second Stage results (Factor Risk Premia and T-stats)
SecondStage <- data.frame(Lambda = Lambda,
                          Tstat = Tstat,
                          Tstat_HAC = Tstat_NW,
                          Tstat_Shanken = Tstat_Shanken)
# Transpose the data frame to have statistics as rows and factors as columns
SecondStage <- t(SecondStage)
# Set row names to statistic names
rownames(SecondStage) <- c('Lambda', 'tstat', 't-stat HAC', 't-stat Shanken')
# Set column names to factor names
colnames(SecondStage) <- c('Mkt', 'C', 'T') # Corrected column names to match factors

# --- Write Results to CSV ---
# Write the First Stage results table to a CSV file
write.csv(FirstStageReg, 'FirstStage_T.csv', row.names = TRUE) # Corrected filename
# Write the Second Stage results table to a CSV file
write.csv(SecondStage, 'SecondStage_T.csv', row.names = TRUE) # Corrected filename


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
