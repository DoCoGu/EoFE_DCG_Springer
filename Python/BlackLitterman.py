# BlackLitterman.py
# This script implements the Black-Litterman model for portfolio allocation,
# blending market equilibrium returns with investor views in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
# Assuming getuncef.py is available in the same directory or Python path
# Make sure getuncef.py is in the same directory as this script, or in your Python path
from getuncef import getuncef

# --- Data Loading and Preparation ---
# Load asset returns from an Excel file (assuming 5 assets)
# The data file is expected to be in the '../Data/' folder relative to this script
ret = pd.read_excel('../Data/Returns.xlsx')
# Load Fama-French factors from an Excel file
# The data file is expected to be in the '../Data/' folder relative to this script
ff = pd.read_excel('../Data/FF_Factors.xlsx')

# Add Market Excess Return and Risk-Free Rate to the returns DataFrame
# Select 'Mkt-RF' and 'RF' columns from the second row onwards (index 1) of the FF DataFrame
# Use .values to get numpy arrays and divide by 100 to convert from percentage to decimal
ret['Mkt-RF'] = ff['Mkt-RF'][1:].values / 100
ret['rf'] = ff['RF'][1:].values / 100

# Calculate excess returns for assets (asset return - risk-free rate)
# Select asset return columns from the 'ret' DataFrame (all except the first, which is Date)
# Exclude the newly added 'Mkt-RF' and 'rf' columns as well
# Use .values to get numpy array and subtract the risk-free rate (reshaped to match dimensions)
exret = ret.iloc[:, 1:-2].values - ret['rf'].values[:, None]  # Excess returns
# Extract Market Excess Return
exmkt = ret['Mkt-RF'].values  # Excess market returns

# --- Compute the Implied Equilibrium Returns (Pi) using CAPM ---
# This step estimates CAPM beta for each asset to derive equilibrium returns based on the market excess return.

T, n = exret.shape # T = number of observations, n = number of assets
# Initialize matrix to store CAPM betas (intercept and slope)
beta = np.empty((2, n))

# Estimate CAPM alpha and beta for each asset by regressing excess asset returns on market excess return
for i in range(n):
    # Add a constant for the intercept (alpha) to the independent variable (exmkt)
    model = sm.OLS(exret[:, i], sm.add_constant(exmkt))
    # Fit the OLS model
    results = model.fit()
    # Store alpha and beta coefficients (results.params[0] is alpha, results.params[1] is beta)
    beta[:, i] = results.params

# Calculate the mean excess return of assets
# Calculate the mean along axis 0 (columns) to get a vector of mean excess returns for each asset
muExret = exret.mean(axis=0)

# Calculate the covariance matrix of the excess returns
# rowvar=False computes the covariance between columns (assets)
Sigma = np.cov(exret, rowvar=False)  # Covariance matrix of the excess returns

# Calculate the implied market equilibrium excess returns (Pi)
# CORRECTED: This calculation matches the logic in your original BlackLitterman.py file:
# Pi = beta * mean_asset_excess_return (element-wise)
# Use the beta coefficient (second row of the beta matrix, index 1 in Python)
Pi = beta[1, :] * muExret  # Market equilibrium returns
# Reshape Pi to be a column vector for matrix calculations later
Pi = Pi[:, None]  # Reshape for calculations

# Uncertainty factor (tau) - often set to a small value or 1/T as defined by Litterman and He (1999)
tau = 1 / T

# --- Define the Investor Views (Q and P) ---
# Define the investor's specific views on expected returns (Q) and the assets involved in those views (P).

# Q is the vector of expected returns on the views, kx1 (k is the number of views)
# Example views:
Q = np.array([0.04 / 12,    # View 1: Expected return for a combination of assets (monthly)
              0.02 / 12,    # View 2: Expected return for another combination (monthly)
              0.10 / 12])[:, None] # View 3: Expected return for a single asset (monthly)

# P is the pick matrix, kxn (k views, n assets)
# Defines which assets are involved in each view and their weights
# Example P matrix corresponding to the example Q views:
# View 1: Asset 1 - Asset 3
# View 2: Asset 4 - Asset 2
# View 3: Asset 5
P = np.array([[1, 0, -1, 0, 0],
              [0, -1, 0, 1, 0],
              [0, 0, 0, 0, 1]])

# Omega is the covariance matrix of the views, kxk
# Represents the uncertainty in the investor's views. It's typically diagonal
# with variances of the view errors. A common approach is P @ (tau * Sigma) @ P.T
Omega = P @ (tau * Sigma) @ P.T

# --- Blend the Equilibrium returns (Pi) with the Views (Q) ---
# Combine the market-implied equilibrium returns with the investor's views using the Black-Litterman formula.

# Calculate the Black-Litterman covariance matrix (SigmaBL)
# Formula: [(tau*Sigma)^-1 + P'*Omega^-1*P]^-1
SigmaBL = np.linalg.inv(np.linalg.inv(tau * Sigma) + P.T @ np.linalg.inv(Omega) @ P)

# Calculate the Black-Litterman expected returns (muBL)
# Formula: SigmaBL * [(tau*Sigma)^-1*Pi + P'*Omega^-1*Q]
muBL = SigmaBL @ (np.linalg.inv(tau * Sigma) @ Pi + P.T @ np.linalg.inv(Omega) @ Q)

# --- Display and Save Expected Returns Comparison ---
# Create a table to compare historical mean excess returns, implied equilibrium returns (Pi), and Black-Litterman blended returns (muBL), and save it to a CSV file.

# Create a pandas DataFrame to display the expected returns comparison
# Concatenate the mean historical excess returns (reshaped to column vector), Pi, and muBL horizontally
# Round the values to 4 decimal places
# Use the original asset names from the 'ret' DataFrame (columns 1 to -2, excluding Date, Mkt-RF, rf) as index
muTable = pd.DataFrame(
    np.round(np.hstack([muExret[:, None], Pi, muBL]), 4),
    columns=['exret', 'Pi', 'muBL'],
    index=ret.columns[1:-2] # Use the original asset names as index
)
print("Comparison of Expected Returns:")
print(muTable)

# Write the table to a CSV file
# The file will be saved in the current working directory (where this script is run)
# index_label='Asset' sets the name for the index column in the CSV
muTable.to_csv('muTable.csv', index_label='Asset')

# --- Implement the two Efficient Frontiers (Mean-Variance and Black-Litterman) ---
# Calculate the efficient frontiers for both the standard Mean-Variance approach (using historical means) and the Black-Litterman approach (using blended returns). This step requires the `getuncef` helper function.

# Define a range of target expected returns for plotting the efficient frontiers
muR = np.arange(0, 0.015, 0.00001)

# Calculate the Mean-Variance efficient frontier using the mean historical excess returns (muExret) and covariance (Sigma)
# getuncef is a helper function assumed to be available (make sure getuncef.py is accessible)
# Pass muExret as a column vector
OptSigma, w = getuncef(muExret[:,None], Sigma, muR)

# Calculate the Black-Litterman efficient frontier using the blended returns (muBL) and adjusted covariance (Sigma + SigmaBL)
# The adjusted covariance matrix is Sigma + SigmaBL (SigmaBL represents the uncertainty in the blended returns)
OptSigmaBL, wBL = getuncef(muBL, Sigma+SigmaBL, muR)

# Create a pandas DataFrame to store the efficient frontier data (returns and standard deviations)
# Flatten OptSigma and OptSigmaBL to ensure they are 1D arrays for DataFrame creation
ef_table = pd.DataFrame(
    {
        'r': muR.flatten(), # Target returns
        'sigma': np.sqrt(OptSigma.flatten()), # Standard deviation for M-V frontier (sqrt of variance)
        'sigmaBL': np.sqrt(OptSigmaBL.flatten()) # Standard deviation for B-L frontier (sqrt of variance)
    }
)

# --- Plot the Efficient Frontiers ---
# Visualize the Mean-Variance and Black-Litterman efficient frontiers.

# Create a figure for the plot
plt.figure(figsize=(10, 6))

# Plot the Mean-Variance efficient frontier
plt.plot(ef_table['sigma'], ef_table['r'], '-b', label='M-V', linewidth=1.5)
# Plot the Black-Litterman efficient frontier
plt.plot(ef_table['sigmaBL'], ef_table['r'], '-r', label='B-L', linewidth=1.5)

# Configure plot appearance
plt.title('Efficient Frontier', fontsize=16)
plt.xlim([0, 0.1]) # Set x-axis limits for better visualization
plt.xlabel('Portfolio Risk', fontsize=16)
plt.ylabel('Portfolio Expected Return', fontsize=16)
plt.legend(loc='upper left', fontsize=14) # Add legend to identify the frontiers
plt.grid(True) # Add a grid for better readability

# Save the plot as a PNG file
# The file will be saved in the current working directory (where this script is run)
plt.savefig("BL.png",format='png')

# plt.show() # Uncomment this line to display the plot directly when running the script

# --- Optimal weights for a specific target return ---
# Calculate the optimal portfolio weights for both the Mean-Variance and Black-Litterman models at a specified target expected return.

# Define a target expected return (e.g., 0.0060)
target_return = 0.0060

# Find the index corresponding to the target return in the efficient frontier table
# Use np.isclose for robust floating point comparison
target_return_idx = np.where(np.isclose(ef_table['r'], target_return))[0]

# Check if the target return was found
if len(target_return_idx) > 0:
    # Extract the weights for the target return from the calculated weight matrices
    # w and wBL are expected to have shape (number of assets, number of target returns)
    # Use .flatten() to ensure the weights are 1D arrays if they were not already
    Weights_mv = w[:, target_return_idx].flatten() # M-V weights
    Weights_bl = wBL[:, target_return_idx].flatten() # B-L weights

    # Create a pandas DataFrame to store the optimal weights for both models
    # Use np.vstack to stack the weight vectors as rows
    # Use the original asset names from the 'ret' DataFrame (columns 1 to 5) as columns
    # Define row names for the two models
    Weights = pd.DataFrame(
        np.vstack([Weights_mv, Weights_bl]),
        columns=ret.columns[1:6], # Use columns 1 to 5 for asset names
        index=["M-V", "B-L"]
    )
    print(f"\nOptimal Weights for Target Return {target_return*100:.2f}%:")
    print(Weights)

    # Write the weights table to an Excel file
    # The file will be saved in the current working directory (where this script is run)
    # sheet_name="Weights" specifies the sheet name in the Excel file
    Weights.to_excel("Weights_BL.xlsx", sheet_name="Weights")
else:
    print(f"\nTarget return {target_return*100:.2f}% not found in the calculated efficient frontier range.")

# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
