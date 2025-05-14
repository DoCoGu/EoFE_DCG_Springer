# CAPM.py
# This script tests the single-factor Capital Asset Pricing Model (CAPM)
# for a set of assets by regressing excess asset returns on the market
# excess return in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import statsmodels.api as sm

# --- Data Loading and Preparation ---
# Import asset returns from an Excel file
Ret = pd.read_excel('../Data/Returns.xlsx')
# Import Fama-French Factors from an Excel file
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

# Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
# Convert from percentage to decimal, starting from the second row (index 1)
Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Reshape Rf to match the shape of asset returns for subtraction
ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# --- Testing the CAPM (One-Factor Model) ---
# Get the size of the excess returns data
nS, mS = ExRet.shape  # nS = number of observations, mS = number of assets

# Initialize arrays to store regression results
# Each array stores the coefficient and its t-statistic for each asset
alpha = np.full((2, mS), np.nan)  # Alpha coefficient and its t-statistic
beta = np.full((2, mS), np.nan)   # Beta coefficient and its t-statistic
R2 = np.full((1, mS), np.nan)    # Adjusted R-squared

# Loop through each asset to perform the CAPM regression
for i in range(mS):
    # Define independent variable (Market Excess Return) and add a constant for the intercept (alpha)
    X = sm.add_constant(Mkt)
    # Define dependent variable (Excess Return for the current asset)
    y = ExRet[:, i]

    # Fit linear model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
    model = sm.OLS(y, X)
    results = model.fit()

    # Extract and store the regression results (coefficients and t-statistics)
    alpha[0, i] = round(results.params[0], 5)  # Alpha coefficient
    alpha[1, i] = round(results.tvalues[0], 5)  # Alpha t-statistic
    beta[0, i] = round(results.params[1], 5)  # Beta coefficient
    beta[1, i] = round(results.tvalues[1], 5)  # Beta t-statistic
    R2[0, i] = round(results.rsquared_adj, 5)  # Adjusted R-squared

# --- Display and Save Results ---
# Create a pandas DataFrame to display the regression results
Stocks = pd.DataFrame(
    data=np.concatenate((alpha, beta, R2)),  # Concatenate results vertically
    columns=Ret.iloc[:, 1:].columns,  # Use original asset names as columns
    index=['alpha', '(alpha t-stat)', 'beta', ' (beta t-stat)', 'Adj. R2']  # Row names
)
print("CAPM Results:")
print(Stocks)

# Write the results table to a CSV file
Stocks.to_csv('CAPM_Stock.csv', index_label='Metric')


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
