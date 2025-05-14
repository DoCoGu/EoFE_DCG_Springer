# MV.py
# This script implements the Mean-Variance portfolio optimization for a set
# of assets and plots the efficient frontier in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Clear the workspace (close all existing plot figures)
plt.close('all')

# --- Data Loading and Preparation ---
# Read returns data from an Excel file
# Using engine='openpyxl' is recommended for .xlsx files
Ret = pd.read_excel('../Data/Returns.xlsx', engine='openpyxl', sheet_name=0)

# Extract return data (assuming it starts from the 2nd column, index 1)
R = Ret.iloc[:, 1:].values
N = R.shape[1] # Number of assets

# Calculate mean returns (expected returns)
z = np.mean(R, axis=0).reshape(-1, 1) # Reshape to a column vector
# Calculate standard deviations (risk)
sig = np.std(R, axis=0)
# Calculate the covariance matrix of returns (rowvar=False assumes columns are variables)
V = np.cov(R, rowvar=False)
# Calculate the inverse of the covariance matrix
V1 = np.linalg.inv(V)

# --- Calculate Parameters for the Efficient Frontier ---
# Calculate the parameters A, B, C, and D from the portfolio theory
# A = z' * V1 * z
A = z.T @ V1 @ z
# B = 1' * V1 * z
B = np.ones((N, 1)).T @ V1 @ z
# C = 1' * V1 * 1
C = np.ones((1, N)) @ V1 @ np.ones((N, 1))
# D = A*C - B**2
D = A * C - B**2

# Define a range of target portfolio expected returns for plotting the efficient frontier
# Using arange for a specific step size
mu_p = np.arange(0, 0.0151, 0.0001)

# Calculate the variance of the efficient frontier for the range of expected returns
# Formula: (1/D) * (C*mu_p**2 - 2*B*mu_p + A)
sig2_p = (1 / D) * (C * mu_p**2 - 2 * B * mu_p + A)
# Calculate the standard deviation (risk) of the efficient frontier
sig_p = np.sqrt(sig2_p)

# --- Plotting the Efficient Frontier ---
# Create a figure
p = plt.figure(1)
# Maximize the figure window (this might behave differently depending on the environment)
p.canvas.manager.full_screen_toggle() # Maximizes the figure window

# Plot the efficient frontier (Risk vs. Expected Return)
plt.plot(sig_p.transpose(), mu_p, '-k', linewidth=1.5) # Plot the efficient frontier line
plt.axhline(y=B/C, linestyle='--', color='k', linewidth=1.5) # Plot horizontal line at GMVP expected return
plt.axvline(1/np.sqrt(C), linestyle='--', color='k', linewidth=1.5) # Plot vertical line at GMVP standard deviation

# Plot scatter points for individual stocks and the Global Minimum Variance Portfolio (GMVP)
plt.scatter(sig, z.flatten(), color='k', marker='o', label='Stocks', alpha=1) # Scatter plot for stocks
plt.scatter(1/np.sqrt(C), B/C, color='k', marker='o', label='MVP', alpha=1) # Scatter plot for GMVP

# Set title, limits, and labels
plt.title('Efficient Frontier', fontsize=16)
plt.xlim(0, 0.1) # Set the x-axis limits
plt.xlabel('Portfolio Risk', fontsize=16)
plt.ylabel('Portfolio Expected Return', fontsize=16)

# Add legend
# Use the specific legend entries and location from the original script
plt.legend(["Efficient frontier", "B/C", "1/sqrt(C)", "Stocks", "MVP"], loc='upper left', fontsize=14)

# Save the figure as an EPS file
plt.savefig('MV.eps', format='eps')

# plt.show() # Uncomment to display the plot

# --- Optimal weights Calculation ---
# Calculate parameters for portfolio weights formula on the efficient frontier
# g and h are vectors that define the linear relationship between weights and expected return
g = 1 / D * (A * (V1 @ np.ones((N, 1))) - B * (V1 @ z))
h = 1 / D * (C * (V1 @ z) - B * (V1 @ np.ones((N, 1))))

# Define a range of target expected returns for which to calculate optimal weights
mup = np.linspace(0.001, 0.05, 10) # Objective return of the portfolio (i.e. what we would like to obtain)
# Calculate optimal weights for the target expected returns
# Formula: w_p = g + h * mu_p
wp = g + h * mup

# --- Global minimum variance portfolio weights ---
# Calculate the weights for the Global Minimum Variance Portfolio (GMVP)
# The GMVP has an expected return of B/C
w_mvp = g + h * (B/C)

# --- Display and Save Weights ---
# Create a pandas DataFrame to store the calculated weights
# Columns are GMVP and target return percentages, index is asset names
columns = ["GMVP"] + [str(round(mu * 100, 2)) + "%" for mu in mup] # Column names
index = Ret.columns[1:] # Row names (asset names)
Weights = pd.DataFrame(index=index, columns=columns, dtype=float) # Initialize empty DataFrame

# Assign values to the table
Weights["GMVP"] = w_mvp.flatten() # Assign GMVP weights
Weights[columns[1:]] = wp # Assign optimal weights for target returns

# Write weights to Excel file
# Save to "Weights_MV.xlsx" without specifying sheet_name (uses default) and with index=True (default)
Weights.to_excel("Weights_MV.xlsx")

# --- Plotting Portfolio Allocation ---
weight_bar = plt.figure(2)
bar_width = 0.8 # Define a width for the bars
bottom = np.zeros(len(mup)) # Initialize the bottom position for stacking

# Iterate through each asset (row of wp) to plot stacked bars
for i in range(N):
    # Plot the bar segment for the current asset's weights across all target returns
    # wp has shape (N, len(mup)), so wp[i, :] gives the weights for asset i across all target returns
    plt.bar(np.arange(len(mup)), wp[i, :], bar_width, bottom=bottom, label=Ret.columns[i+1])
    # Update the bottom position for the next asset's stack
    bottom += wp[i, :]

# Configure plot appearance
# Corrected legend location from 'northwest' to 'upper left' (a valid option)
plt.legend(loc='upper left') # Add legend with asset names
plt.xticks(np.arange(len(mup)), [f"{round(ret*100, 2)}%" for ret in mup]) # Set x-axis ticks and labels to target returns
plt.ylabel("Weight (%)")
plt.xlabel("Expected return (%)")
plt.title("Portfolio allocation")
# Set figure position and size (this might behave differently depending on the environment)
# weight_bar.set_size_inches(8, 4) # Example of setting size

plt.tight_layout() # Adjust layout to prevent labels overlapping
plt.savefig("MV_p.eps", format='eps') # Save the plot as an EPS file
# plt.show() # Uncomment to display the plot


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
