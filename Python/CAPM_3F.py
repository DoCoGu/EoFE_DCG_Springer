# CAPM_3F.py
# This script tests the Fama-French three-factor model for a set of assets
# by regressing excess asset returns on the market, SMB, and HML factors in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import statsmodels.api as sm

# --- Data Loading and Preparation ---
# Load asset returns from an Excel file
Ret = pd.read_excel('../Data/Returns.xlsx')
# Load Fama-French Factors from an Excel file
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

# Extract Fama-French factors (Mkt-RF, SMB, HML)
# Convert from percentage to decimal, starting from the second row (index 1)
FF = np.array(Factors.iloc[1:, 1:4]) / 100

# Extract Risk-Free Rate (RF)
# Convert from percentage to decimal, starting from the second row (index 1)
Rf = np.array(Factors.iloc[1:, 4]) / 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Reshape Rf to match the shape of asset returns for subtraction
ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# --- Testing the Fama-French 3-Factor Model ---
# Get the size of the excess returns data
nS, mS = ExRet.shape  # nS = number of observations, mS = number of assets

# Initialize arrays to store regression results for the 3-Factor model
# Each array stores the coefficient and its t-statistic for each asset
alpha = np.empty((2, mS))  # Alpha coefficient and its t-statistic
beta1 = np.empty((2, mS))  # Beta for Market factor and its t-statistic
beta2 = np.empty((2, mS))  # Beta for SMB factor and its t-statistic
beta3 = np.empty((2, mS))  # Beta for HML factor and its t-statistic
R2 = np.empty((1, mS))   # Adjusted R-squared

# Loop through each asset to perform the Fama-French 3-Factor regression
for i in range(mS):
    # Define independent variables (FF factors) and add a constant for the intercept (alpha)
    X = sm.add_constant(FF)
    # Define dependent variable (Excess Return for the current asset)
    y = ExRet[:, i]

    # Fit linear model: ExRet_i = alpha_i + beta1_i*Mkt + beta2_i*SMB + beta3_i*HML + epsilon_i
    model = sm.OLS(y, X)
    results = model.fit()

    # Extract and store the regression results (coefficients and t-statistics)
    alpha[0, i] = round(results.params[0], 5)  # Alpha coefficient
    alpha[1, i] = round(results.tvalues[0], 5)  # Alpha t-statistic
    beta1[0, i] = round(results.params[1], 5)  # Market Beta coefficient
    beta1[1, i] = round(results.tvalues[1], 5)  # Market Beta t-statistic
    beta2[0, i] = round(results.params[2], 5)  # SMB Beta coefficient
    beta2[1, i] = round(results.tvalues[2], 5)  # SMB Beta t-statistic
    beta3[0, i] = round(results.params[3], 5)  # HML Beta coefficient
    beta3[1, i] = round(results.tvalues[3], 5)  # HML Beta t-statistic
    R2[0, i] = round(results.rsquared_adj, 5)  # Adjusted R-squared

# --- Display and Save Results ---
# Create a pandas DataFrame to display the regression results
Stocks = pd.DataFrame(
    data=np.concatenate((alpha, beta1, beta2, beta3, R2)),  # Concatenate results vertically
    columns=Ret.iloc[:, 1:].columns,  # Use original asset names as columns
    index=['alpha', '(alpha t-stat)', 'beta_mkt', 'beta_mkt (t-stat)',
           'beta_smb', '(beta_smb t-stat)', 'beta_hml', '(beta_hml t-stat)', 'Adj. R2']  # Row names
)
print("Fama-French 3-Factor Model Results:")
print(Stocks)

# Write the results table to a CSV file
Stocks.to_csv('CAPM_3F_Stock.csv', index_label='Metric')


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
