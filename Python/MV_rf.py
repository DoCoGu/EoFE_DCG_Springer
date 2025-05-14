# MV_rf.py
# This script implements Mean-Variance portfolio optimization with a
# risk-free asset, deriving and plotting the Capital Market Line (CML) in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Assuming getuncef.py is available for portfolio calculations if needed
# from getuncef import getuncef # Uncomment if getuncef is used

# --- Data Loading and Preparation ---
# Import data from Excel files
# Loading data assuming files are in the '../Data/' directory relative to this script
Ret = pd.read_excel('../Data/Returns.xlsx')
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

# Extract risk-free rate (assuming it's in the 5th column, index 4, starting from the 2nd row, index 1)
# Convert from annual percentage to monthly decimal
Rf = np.array(Factors.iloc[1:, 4]) / 1200

# Extract asset return data (assuming it starts from the 2nd column, index 1)
R = np.array(Ret.iloc[:, 1:])
N = R.shape[1] # Number of assets

# Calculate mean returns (expected returns)
# Using reshape(-1, 1) to ensure it's a column vector for robust matrix operations
z = np.mean(R, axis=0).reshape(-1, 1) # Reshape to a column vector
# Calculate standard deviations (risk)
sig = np.std(R, axis=0)
# Calculate the covariance matrix of returns (rowvar=False assumes columns are variables)
V = np.cov(R, rowvar=False)
# Calculate the inverse of the covariance matrix
V1 = np.linalg.inv(V)

# --- Calculate Parameters for the Capital Market Line (CML) ---
# Calculate H, which is related to the squared Sharpe Ratio of the tangency portfolio
# H = (z - rf*1)' * inv(S) * (z - rf*1)
# Using proper matrix multiplication with z as a column vector
H = (z - np.mean(Rf)).T @ V1 @ (z - np.mean(Rf))
sqrtH = np.sqrt(H) # Related to the slope of the CML

# Define a range of target portfolio expected returns for plotting the CML
mu_p = np.linspace(np.mean(Rf), 0.009, 100)

# Calculate the variance of portfolios on the Capital Market Line
# Formula: (1/H) * (mu_p - rf)^2
sig2_p = (1 / H) * (mu_p - np.mean(Rf))**2
# Calculate the standard deviation (risk) of portfolios on the CML
sig_p = np.sqrt(sig2_p)

# Calculate the expected return and standard deviation of the Tangency Portfolio
# mu_t = rf + H / (1' * inv(S) * (z - rf*1))
# Using proper matrix operations
mu_t = np.mean(Rf) + H / (np.ones((1, N)) @ V1 @ (z - np.ones((N, 1)) * np.mean(Rf)))
# sig_t = sqrt(H) / (1' * inv(S) * (z - rf*1))
# Using proper matrix operations
sig_t = sqrtH / (np.ones((1, N)) @ V1 @ (z - np.ones((N, 1)) * np.mean(Rf)))


# --- Calculate Parameters for the Risky Asset Efficient Frontier (for comparison) ---
# Recalculate parameters A, B, C, D for the risky asset efficient frontier
# Using proper matrix operations with z as a column vector
A = z.T @ V1 @ z
B = np.ones((1, N)) @ V1 @ z
C = np.ones((1, N)) @ V1 @ np.ones((N, 1))
D = A * C - B**2

# Define a range of target expected returns for the risky asset efficient frontier
mu_p_risky = np.linspace(0.001, 0.008, 100)

# Calculate the variance of the risky asset efficient frontier
sig2_p_risky = (1 / D) * (C * mu_p_risky**2 - 2 * B * mu_p_risky + A)
sig_pp_risky = np.sqrt(sig2_p_risky)
# Set standard deviations below a certain threshold to NaN for cleaner plotting
# Using np.where for robust conditional assignment
sig_pp_risky = np.where(mu_p_risky < 0.005, np.nan, sig_pp_risky)
# x = np.vstack((sig_pp_risky, mu_p_risky)) # Combine for potential use (though not directly used in plot)


# --- Plotting the Capital Market Line and Efficient Frontier ---
fig = plt.figure(2) # Use figure 2 as in the original script
# Plot the Capital Market Line (CML)
plt.plot(sig_p.flatten(), np.mean(Rf) + sqrtH.flatten() * sig_p.flatten(), '-b', linewidth=1.5) # Ensure 1D for plotting
# Plot the risky asset efficient frontier
plt.plot(sig_pp_risky.flatten(), mu_p_risky.flatten(), '-k', linewidth=1.5) # Ensure 1D for plotting
# Scatter plot of the Tangency Portfolio
plt.scatter(sig_t.flatten(), mu_t.flatten(), marker='o', color='k') # Ensure 1D for plotting

# Configure plot appearance
plt.title('Capital market line', fontsize=16) # Title as in original
plt.xlabel('Portfolio Risk', fontsize=16) # X-label as in original
plt.ylabel('Portfolio Expected Return', fontsize=16) # Y-label as in original
# Use the exact legend entries from the original script
# Note: Some of these labels ('Portfolio', 'Stocks', 'Autoupdate', 'off') might not directly correspond
# to plotted elements in this figure, but are included to match the original legend call.
plt.legend(["Efficient frontier", "Portfolio", "Tangency port.", "Stocks", "Autoupdate", "off"], fontsize=16, loc="upper left")
plt.xlim([0, 0.05]) # Set x-axis limits as in original
plt.ylim([np.mean(Rf), 0.009]) # Set y-axis limits, starting from the risk-free rate

# Add dashed lines to highlight the Tangency Portfolio's position
plt.axhline(mu_t.flatten(), linestyle='--', color='k', linewidth=1.5) # Ensure 1D for plotting
plt.axvline(sig_t.flatten(), linestyle='--', color='k', linewidth=1.5) # Ensure 1D for plotting

# Save the plot as an EPS file
plt.savefig("MV_rf.eps", format="eps")
# plt.show() # Uncomment to display the plot

# --- Optimal weights Calculation ---
# Define a range of target expected returns for which to calculate optimal weights
mup_target = np.linspace(0.001, 0.05, 10)
# Calculate optimal weights for portfolios on the CML (combination of risk-free and tangency)
# Formula: w_p = inv(S) * (z - rf*1) * (mu_p - rf) / H
# Using proper matrix operations
wp = V1 @ (z - np.ones((N, 1)) * np.mean(Rf)) @ ((mup_target - np.mean(Rf)) / H)

# Calculate the weights for the Tangency Portfolio (risky assets only)
# Formula: w_T = inv(S) * (z - rf*1) / (1' * inv(S) * (z - rf*1))
# Using proper matrix operations
wT = V1 @ (z - np.ones((N, 1)) * np.mean(Rf)) / (np.ones((1, N)) @ V1 @ (z - np.ones((N, 1)) * np.mean(Rf)))

# Calculate the weight allocated to the risk-free asset for each portfolio (Tangency and CML)
# rf_weight = 1 - sum of risky asset weights in each portfolio
# Sum of risky weights for Tangency Portfolio (wT)
sum_wT = np.sum(wT)
# Sum of risky weights for CML portfolios (wp) along axis 0
sum_wp = np.sum(wp, axis=0)
# Combine into a single row vector: [1 - sum(wT), 1 - sum(wp[:, 0]), ..., 1 - sum(wp[:, 9])]
rf_weight_row = np.hstack([1 - sum_wT, 1 - sum_wp])
# Reshape to a row vector of shape (1, 1 + len(mup_target))
rf_weight_row = rf_weight_row.reshape(1, -1)

# --- Display and Save Weights ---
# Create the full weight data matrix by vertically stacking risky weights and the risk-free weight row
# Risky weights matrix: horizontally stack wT (N, 1) and wp (N, len(mup_target)) -> shape (N, 1 + len(mup_target))
risky_weights_matrix = np.hstack([wT, wp])
# Vertically stack the risky weights matrix and the risk-free weight row
weight_data = np.vstack([risky_weights_matrix, rf_weight_row])

# Round the weight data to 3 decimal places as in the original script
weight_data = np.round(weight_data, 3)

# Define column names (Tangency Portfolio + target return percentages)
columns = ["wT"] + [f"{round(mup_val * 100, 2)}%" for mup_val in mup_target]
# Define row names (asset names + Rf)
index = Ret.columns[1:].tolist() + ["Rf"]

# Create a pandas DataFrame to store the calculated weights
Weights = pd.DataFrame(weight_data, columns=columns, index=index)

print("Optimal Weights (Tangency and CML):")
print(Weights)

# Write the weights table to an Excel file
# Save to "Weights_MVrf.xlsx" with index=True and header=True (default behavior)
Weights.to_excel("Weights_MVrf.xlsx", index=True, header=True)


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
