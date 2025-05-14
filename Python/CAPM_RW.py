# CAPM_RW.py
# This script implements a rolling-window regression to estimate time-varying
# CAPM alpha and beta for a set of assets using Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

# --- Data Loading and Preparation ---
# Load asset returns from an Excel file
Ret = pd.read_excel('../Data/Returns.xlsx')
# Load Fama-French Factors from an Excel file
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

# Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
# Convert from percentage to decimal, starting from the second row (index 1)
Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Reshape Rf to match the shape of asset returns for subtraction
ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# --- Rolling Regression Setup ---
# Define the size of the rolling window
rw = 40

# Get the size of the excess returns data
nS, mS = ExRet.shape  # nS = total number of observations, mS = number of assets

# Initialize arrays to store rolling regression results
# alpha and beta will store the coefficients for each window and asset
alpha = np.empty((nS - rw + 1, mS))
beta = np.empty((nS - rw + 1, mS))

# --- Perform Rolling Window Regression ---
# Outer loop iterates through the starting point of each window
for t in range(nS - rw + 1):
    # Inner loop iterates through each asset
    for i in range(mS):
        # Define the current window for the independent variable (Market Excess Return)
        X_window = Mkt[t:t + rw - 1]
        # Define the current window for the dependent variable (Asset Excess Return)
        y_window = ExRet[t:t + rw - 1, i]

        # Add a constant to the independent variable for the intercept (alpha) in OLS
        X_window = sm.add_constant(X_window)

        # Fit the linear regression model for the current window and asset
        # Model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
        model = sm.OLS(y_window, X_window)
        results = model.fit()

        # Extract and store the regression results (alpha and beta) for the current window
        alpha[t, i] = round(results.params[0], 5)  # Alpha coefficient
        beta[t, i] = round(results.params[1], 5)  # Beta coefficient

# --- Plotting Rolling Estimates ---
# Prepare dates for plotting the rolling estimates
# Create a sequence of dates corresponding to the end of each rolling window
# Use the date column from the original Returns DataFrame
Dates_for_plot = pd.to_datetime(Ret.iloc[rw - 1:, 0]).tolist()
# Select 20 equally spaced dates for x-axis ticks
Dates_for_plot = Dates_for_plot[::len(Dates_for_plot)//20 + 1]
# Format the dates as 'Month-Year' strings
Dates_for_plot = [d.strftime('%b-%Y') for d in Dates_for_plot]

# alpha plot
alpha_plot = plt.figure(1)
# Plot alpha for each asset over time
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 0], 'r', linewidth=1.5, label=Ret.columns[1])
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 1], 'k', linewidth=1.5, label=Ret.columns[2])
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 2], 'b', linewidth=1.5, label=Ret.columns[3])
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 3], 'g', linewidth=1.5, label=Ret.columns[4])
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 4], 'm', linewidth=1.5, label=Ret.columns[5])
# Add a horizontal line at zero for reference
plt.plot(np.arange(1, nS - rw + 2), np.zeros(nS - rw + 1), '--k')

# Configure plot appearance
plt.xticks(np.linspace(1, nS - rw + 1, len(Dates_for_plot)), Dates_for_plot, rotation=45, fontsize=14) # Set x-axis ticks and labels
plt.xlim([1, nS - rw + 1]) # Set x-axis limits
plt.legend(loc='best') # Add legend
plt.title('Rolling Alpha (CAPM)')
plt.xlabel('Time')
plt.ylabel('Alpha')
plt.tight_layout() # Adjust layout to prevent labels overlapping
plt.savefig('alpha_S.png') # Save the plot as a PNG file
# plt.show() # Uncomment to display the plot

# beta plot
beta_plot = plt.figure(2)
# Plot beta for each asset over time
plt.plot(np.arange(1, nS - rw + 2), beta[:, 0], 'r', linewidth=1.5, label=Ret.columns[1])
plt.plot(np.arange(1, nS - rw + 2), beta[:, 1], 'k', linewidth=1.5, label=Ret.columns[2])
plt.plot(np.arange(1, nS - rw + 2), beta[:, 2], 'b', linewidth=1.5, label=Ret.columns[3])
plt.plot(np.arange(1, nS - rw + 2), beta[:, 3], 'g', linewidth=1.5, label=Ret.columns[4])
plt.plot(np.arange(1, nS - rw + 2), beta[:, 4], 'm', linewidth=1.5, label=Ret.columns[5])
# Add a horizontal line at one for reference (CAPM beta is expected to be 1 for the market)
plt.plot(np.arange(1, nS - rw + 2), np.ones(nS - rw + 1), '--k')

# Configure plot appearance
plt.xticks(np.linspace(1, nS - rw + 1, len(Dates_for_plot)), Dates_for_plot, rotation=45, fontsize=14) # Set x-axis ticks and labels
plt.xlim([1, nS - rw + 1]) # Set x-axis limits
plt.legend(loc='best') # Add legend
plt.title('Rolling Beta (CAPM)')
plt.xlabel('Time')
plt.ylabel('Beta')
plt.tight_layout() # Adjust layout to prevent labels overlapping
plt.savefig('beta_S.png') # Save the plot as a PNG file
# plt.show() # Uncomment to display the plot


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
