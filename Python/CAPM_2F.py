# CAPM_2F.py
# This script tests two-factor asset pricing models by adding either VIX
# changes or EPU changes as a second factor to the CAPM in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import statsmodels.api as sm

# --- Data Loading and Preparation ---
# Import Stock Returns from an Excel file
Ret = pd.read_excel('../Data/Returns.xlsx')
# Import Fama-French Factors from an Excel file
Factors = pd.read_excel('../Data/FF_Factors.xlsx')
# Import Uncertainty data (including VIX and EPU) from an Excel file
Uncertainty = pd.read_excel('../Data/Uncertainty.xlsx')

# Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
# Convert from percentage to decimal, starting from the second row (index 1)
Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100

# Calculate Excess Returns for stocks (Stock Return - Risk-Free Rate)
# Reshape Rf to match the shape of stock returns for subtraction
ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# Calculate percentage change in VIX and EPU
# VIX change: (Current VIX / Previous VIX) - 1
VIX = np.array(Uncertainty.iloc[1:, 1].values / Uncertainty.iloc[:-1, 1].values - 1)
# EPU change: (Current EPU / Previous EPU) - 1
EPU = np.array(Uncertainty.iloc[1:, 2].values / Uncertainty.iloc[:-1, 2].values - 1)

# --- Testing the CAPM with VIX (Two-Factor Model) ---
# Get the size of the excess returns data
nS, mS = ExRet.shape # nS = number of observations, mS = number of stocks

# Initialize arrays to store regression results for CAPM + VIX model
# Each array stores the coefficient and its t-statistic for each stock
alpha_VIX = np.empty((2, mS))
beta1_VIX = np.empty((2, mS)) # Beta for Market factor
beta2_VIX = np.empty((2, mS)) # Beta for VIX factor
R2_VIX = np.empty((1, mS)) # Adjusted R-squared

# Loop through each stock to perform the regression with Market and VIX
for i in range(mS):
    # Define independent variables (Market and VIX) and add a constant for the intercept
    X_VIX = sm.add_constant(np.array([Mkt, VIX]).transpose())
    # Define dependent variable (Excess Return for the current stock)
    y = ExRet[:, i]

    # Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * VIX + epsilon_i
    model = sm.OLS(y, X_VIX)
    results = model.fit()

    # Extract and store the regression results (coefficients and t-statistics)
    alpha_VIX[0, i] = round(results.params[0], 5) # Alpha coefficient
    alpha_VIX[1, i] = round(results.tvalues[0], 5) # Alpha t-statistic
    beta1_VIX[0, i] = round(results.params[1], 5) # Market Beta coefficient
    beta1_VIX[1, i] = round(results.tvalues[1], 5) # Market Beta t-statistic
    beta2_VIX[0, i] = round(results.params[2], 5) # VIX Beta coefficient
    beta2_VIX[1, i] = round(results.tvalues[2], 5) # VIX Beta t-statistic
    R2_VIX[0, i] = round(results.rsquared_adj, 5) # Adjusted R-squared

# --- Testing the CAPM with EPU (Two-Factor Model) ---
# Initialize arrays to store regression results for CAPM + EPU model
alpha_EPU = np.empty((2, mS))
beta1_EPU = np.empty((2, mS)) # Beta for Market factor
beta2_EPU = np.empty((2, mS)) # Beta for EPU factor
R2_EPU = np.empty((1, mS)) # Adjusted R-squared

# Loop through each stock to perform the regression with Market and EPU
for i in range(mS):
    # Define independent variables (Market and EPU) and add a constant for the intercept
    X_EPU = sm.add_constant(np.array([Mkt, EPU]).transpose())
    # Define dependent variable (Excess Return for the current stock)
    y = ExRet[:, i]

    # Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * EPU + epsilon_i
    model = sm.OLS(y, X_EPU)
    results = model.fit()

    # Extract and store the regression results (coefficients and t-statistics)
    alpha_EPU[0, i] = round(results.params[0], 5) # Alpha coefficient
    alpha_EPU[1, i] = round(results.tvalues[0], 5) # Alpha t-statistic
    beta1_EPU[0, i] = round(results.params[1], 5) # Market Beta coefficient
    beta1_EPU[1, i] = round(results.tvalues[1], 5) # Market Beta t-statistic
    beta2_EPU[0, i] = round(results.params[2], 5) # EPU Beta coefficient
    beta2_EPU[1, i] = round(results.tvalues[2], 5) # EPU Beta t-statistic
    R2_EPU[0, i] = round(results.rsquared_adj, 5) # Adjusted R-squared

# --- Display and Save Results ---
# Create a pandas DataFrame for the CAPM + VIX results
Stocks_VIX = pd.DataFrame(
    data=np.concatenate((alpha_VIX, beta1_VIX, beta2_VIX, R2_VIX)), # Concatenate results vertically
    columns=Ret.iloc[:, 1:].columns, # Use original stock names as columns
    index=['alpha', '(alpha t-stat)', 'Mkt', ' (Mkt t-stat)', 'VIX', '(VIX t-stat)', 'Adj. R2'] # Row names
)
print("CAPM + VIX Results:")
print(Stocks_VIX)
Stocks_VIX.to_csv('CAPM_Stock_VIX.csv', index_label='Metric') # Write results to CSV

# Create a pandas DataFrame for the CAPM + EPU results
Stocks_EPU = pd.DataFrame(
    data=np.concatenate((alpha_EPU, beta1_EPU, beta2_EPU, R2_EPU)), # Concatenate results vertically
    columns=Ret.iloc[:, 1:].columns, # Use original stock names as columns
    index=['alpha', '(alpha t-stat)', 'Mkt', ' (Mkt t-stat)', 'EPU', '(EPU t-stat)', 'Adj. R2'] # Row names
)
print("\nCAPM + EPU Results:")
print(Stocks_EPU)
Stocks_EPU.to_csv('CAPM_Stock_EPU.csv', index_label='Metric') # Write results to CSV


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
