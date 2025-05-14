# CS_EPU.py
# This script performs Fama-MacBeth cross-sectional regression analysis,
# including Economic Policy Uncertainty (EPU) and consumption growth as factors in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import statsmodels.api as sm
# Assuming getuncef.py is available for portfolio calculations if needed,
# although this script focuses on cross-sectional regressions.
# from getuncef import getuncef # Uncomment if getuncef is used

# --- Data Loading and Preparation ---
# Import Portfolio Returns from an Excel file
Ret = pd.read_excel('../Data/Portfolios.xlsx')
# Import Fama-French Factors from an Excel file
FFactors = pd.read_excel('../Data/FF_Factors.xlsx')
# Import Uncertainty data (including EPU) from an Excel file
Uncertainty = pd.read_excel('../Data/Uncertainty.xlsx')
# Import Consumption data from an Excel file
Consumption = pd.read_excel('../Data/Consumption.xlsx')

# Calculate percentage change in EPU and Consumption
# EPU change: (Current EPU / Previous EPU) - 1
# Using .values to get numpy array and slicing to align data
EPU = Uncertainty.iloc[1:, 2].values / Uncertainty.iloc[:-1, 2].values - 1
# Consumption change: (Current Consumption / Previous Consumption) - 1
C = Consumption.iloc[1:, 1].values / Consumption.iloc[:-1, 1].values - 1

# Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
# Convert from percentage to decimal, starting from the second row (index 1)
Mkt = FFactors.iloc[1:, 1].values / 100
Rf = FFactors.iloc[1:, 4].values / 100

# Calculate Excess Returns for portfolios (Portfolio Return - Risk-Free Rate)
# Note: Using data from row 2 onwards (index 1) for Ret and Rf to match lengths
ExRet = Ret.iloc[1:, 1:].values / 100 - Rf.reshape(-1, 1)

# Combine the factors into a pandas DataFrame
Factors = pd.DataFrame({'Mkt': Mkt, 'C': C, 'EPU': EPU})

# --- First Stage Regression (Time-Series Regressions) ---
# Get the size of the excess returns and factors data
n1, n2 = ExRet.shape # n1 = number of observations, n2 = number of portfolios
nF = Factors.shape[1] # nF = number of factors

# Initialize arrays to store results from the first stage
# CoefAll stores the factor betas for each portfolio (excluding intercept)
CoefAll = np.empty((nF, n2))
# Res stores the residuals from each regression
Res = np.empty((n1, n2))

# Loop through each portfolio to perform time-series regression
for i in range(n2):
    # Regress portfolio excess returns on the factors
    # Model: ExRet_i = alpha_i + beta_Mkt*Mkt + beta_dC*dC + beta_EPU*EPU + epsilon_i
    # Add a constant to the factors for the intercept (alpha) in OLS
    model = sm.OLS(ExRet[:, i], sm.add_constant(Factors)).fit()

    # Store the factor betas (excluding the intercept, which is at index 0)
    CoefAll[:, i] = round(model.params[1:], 5)
    # Store the residuals
    Res[:, i] = round(model.resid, 5)

# Calculate the variance-covariance matrix of the residuals (transpose Res for cov function)
VarCovErr = np.cov(Res.T)

# --- Second-Stage Regression (Cross-Sectional Regressions) ---
# Calculate the average excess return for each portfolio
MeanRet = np.mean(ExRet, axis=0)

# Transpose the factor betas matrix for the cross-sectional regression
Betas = CoefAll.T

# Perform cross-sectional regression of average returns on factor betas
# Model: MeanRet_i = lambda_0 + lambda_Mkt*beta_Mkt_i + lambda_dC*beta_dC_i + lambda_EPU*beta_EPU_i + eta_i
# Using 'Intercept', false assumes lambda_0 = 0 in the cross-sectional model
model = sm.OLS(MeanRet, Betas).fit()

# Extract and store the factor risk premia (Lambdas) and their standard errors
SE = model.bse # Standard Errors of the Lambdas
Lambda = model.params # Factor Risk Premia (Lambdas)
Tstat = Lambda / SE # T-statistics (standard)

# Calculate Shanken-corrected standard errors
Sigma_f = np.cov(Factors, rowvar=False) # Covariance matrix of factors

B = Betas # Rename Betas for clarity in formula
BtB_inv = np.linalg.inv(B.T @ B) # Inverse of (Betas' * Betas)
# Calculate the correction term for Shanken standard errors
correction = 1 + Lambda.T @ np.linalg.inv(Sigma_f) @ Lambda
# Calculate the Shanken-corrected covariance matrix of Lambdas
VarLam = BtB_inv @ B.T @ VarCovErr @ B @ BtB_inv * correction + Sigma_f
# Divide by the number of observations (n1) as per the Fama-MacBeth approach
VarLam = VarLam / n1

SE_Shanken = np.sqrt(np.diag(VarLam)) # Shanken-corrected Standard Errors (diagonal of covariance matrix)
Tstat_Shanken = Lambda / SE_Shanken # Shanken-corrected T-statistics

# --- Time-series of cross-sectional regressions (for Newey-West standard errors) ---
nF = Betas.shape[1] # Number of factors (same as number of betas)
LambdaFull = np.full((n1, nF), np.nan) # Matrix to store Lambdas from each cross-sectional regression

# Loop through each time period (observation)
for j in range(n1):
    # Use returns from a single time period for the cross-sectional regression
    MeanRet_j = ExRet[j, :]
    # Fit cross-sectional model for the current time period (using betas from the first stage)
    model_j = sm.OLS(MeanRet_j, Betas).fit()
    # Store the estimated Lambdas for this time period
    LambdaFull[j, :] = model_j.params

# Calculate the mean of the estimated Lambdas across all time periods
LambdaMean = np.mean(LambdaFull, axis=0)

# Calculate HAC-corrected (Newey-West) standard errors
# This requires the time series of estimated Lambdas (LambdaFull)
X = np.ones((n1, 1)) # Dummy predictor for HAC calculation (intercept only)
hac_cov = np.zeros(nF) # Initialize array for HAC variances

# Loop through each factor to calculate its HAC standard error
for k in range(nF):
    # Use the time series of estimated Lambdas for the current factor
    y = LambdaFull[:, k]
    # Check if y contains non-NaN values before fitting the model
    if not np.isnan(y).all():
        # Fit a dummy regression to use the HAC function from statsmodels
        model_k = sm.OLS(y, X)
        try:
            results_k = model_k.fit(cov_type='HAC', cov_kwds={'maxlags': 2, 'kernel': 'bartlett'})
            hac_cov[k] = results_k.cov_params()[0, 0] # Extract the variance (top-left element)
        except ValueError as e:
             print(f"Could not calculate HAC standard error for factor {k}: {e}")
             hac_cov[k] = np.nan # Set to NaN if calculation fails
    else:
        hac_cov[k] = np.nan # Set to NaN if all Lambdas are NaN

SE_NW = np.sqrt(hac_cov) # Newey-West Standard Errors

# Calculate Newey-West T-statistics (handle division by zero if SE_NW is 0 or NaN)
Tstat_NW = np.divide(LambdaMean, SE_NW, out=np.full_like(LambdaMean, np.nan), where=SE_NW!=0)


# --- Display and Save Results ---
# Get the names of the portfolios from the original data
NamePort = Ret.iloc[:, 1:].columns

# Create a pandas DataFrame for the First Stage results (Factor Betas)
FirstStageReg = pd.DataFrame(
    Betas,
    columns=['Mkt', 'dC', 'EPU'], # Column names for factors
    index=NamePort # Row names for portfolios
)
print("First Stage Regression (Factor Betas):")
print(FirstStageReg)
FirstStageReg.to_csv('FirstStage_EPU.csv', index_label='Portfolio') # Write to CSV

# Create a pandas DataFrame for the Second Stage results (Factor Risk Premia and T-stats)
SecondStage = pd.DataFrame(
    np.vstack([Lambda, Tstat, Tstat_NW, Tstat_Shanken]),
    columns=['Mkt', 'dC', 'EPU'], # Column names for factors
    index=['Lambda', 'tstat', 't-stat HAC', 't-stat Shanken'] # Row names for statistics
)
print("\nSecond Stage Regression (Factor Risk Premia and T-statistics):")
print(SecondStage)
SecondStage.to_csv('SecondStage_EPU.csv', index_label='Statistic') # Write to CSV


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
