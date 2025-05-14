# ES_RegAnalysis.py
# This script performs regression analysis on industry portfolio returns,
# including event studies and calendar effects (weekday and tax days) in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import statsmodels.api as sm

# --- Setting Parameters ---
maxlag = 3  # Maximum lag for lagged returns in the regression
cutoff = 150  # Cutoff value for filtering events (based on column 13 in Events.xlsx)
fD = 3  # Number of future days to create event dummies for
cons = 1 # Flag for including a constant in the regression (1 for include) - Note: statsmodels OLS adds constant separately
America = 1  # Flag to filter events by Zone (1 for America, 0 for All)

# --- Data Loading and Preparation ---
# Import industry portfolio returns from an Excel file
Ret = pd.read_excel("../Data/IndustryPortfolios.xlsx")
# Extract dates from the first column and convert to datetime objects
DatesReturn = pd.to_datetime(Ret.iloc[:, 0], format='%Y%m%d')

# Import events data from an Excel file
Events = pd.read_excel('../Data/Events.xlsx')

# Filter events based on the America flag and casualty cutoff
if America == 0:
    # Filter for all events above the cutoff (assuming casualties are in column 12, event date in column 3)
    Events = Events[Events.iloc[:, 12] > cutoff].iloc[:, 3]
else:
    # Filter for events in Zone 1 (America, assuming Zone is in the last column) above the cutoff
    Events = Events[(Events.iloc[:, 12] > cutoff) & (Events["Zone"] == 1)].iloc[:, 3]

# Extract returns data (assuming it's in the 14th column, index 13)
Ret = Ret.iloc[:, 13].values

# Convert event dates to datetime objects (handling potential NaT from filtering)
Events = pd.to_datetime(Events.values.flatten(), format='%d-%m-%Y', errors='coerce')
Events = Events.dropna() # Remove any dates that failed to convert

# --- Create Dummy Variables ---
# Preallocate array for weekday dummies (Monday to Thursday)
dow = np.full((len(DatesReturn), 4), np.nan)

# Create dummy variables for weekdays (Monday=0, Tuesday=1, ..., Thursday=3 in pandas weekday)
for d in range(4): # Loop for Monday (0) to Thursday (3)
    for i in range(len(DatesReturn)):
        dow[i, d] = 1 if DatesReturn[i].weekday() == d else 0

# Create events dummy variables
# Check if each date in DatesReturn is present in the filtered Events dates
Match = np.isin(DatesReturn, Events)

# Initialize event dummy matrix. Adding extra rows initially to handle shifting.
# Shape: (number of return dates + fD, fD)
EDummy = np.zeros((len(DatesReturn) + fD, fD))

# Create dummies for event days and subsequent days (up to fD)
for i in range(len(DatesReturn)):
    if Match[i]: # If the current date is an event date
        # Set dummies for the next fD days
        for fd_idx in range(fD):
            if (i + 1 + fd_idx) < (len(DatesReturn) + fD): # Ensure index is within bounds
                EDummy[i + 1 + fd_idx, fd_idx] = 1

# Remove the first fD rows as they cannot have event dummies referring to previous events
EDummy = EDummy[fD:, :]

# Create a dummy variable for tax days (Jan 1-5)
T = np.zeros(len(Ret))
dday = DatesReturn.dt.day
dmonth = DatesReturn.dt.month

# Identify dates that are in January and between the 1st and 5th
T[(dmonth == 1) & (dday >= 1) & (dday <= 5)] = 1

# --- Prepare Regression Data ---
# Create lagged returns matrix
# np.roll shifts the array. We roll by lag and set the first 'lag' elements to NaN.
RetLag = np.column_stack([np.roll(Ret, lag) for lag in range(1, maxlag + 1)])
RetLag[:maxlag, :] = np.nan # Set the initial lagged values to NaN

# Combine all independent variables into a single matrix X
# Stack lagged returns, weekday dummies, tax day dummy, and event dummies
X = np.column_stack([RetLag, dow, T, EDummy[:, :fD]])

# Remove rows with NaN values (due to lagged returns) from both X and the dependent variable (Ret)
# This aligns the data for regression
X = X[maxlag:]
Ret_trimmed = Ret[maxlag:]

# Adding a constant for the intercept in statsmodels OLS
X = sm.add_constant(X)

# --- Fit the Linear Regression Model ---
mdl = sm.OLS(Ret_trimmed, X).fit()

# Extract coefficients and standard errors
coeff = np.round(mdl.params, 5)
se = mdl.bse  # Standard errors

# Calculate t-statistics (handle division by zero if se is 0)
Sig = np.divide(coeff, se, out=np.full_like(coeff, np.nan), where=se!=0)

# Combine coefficients and t-statistics for display
res = np.vstack([coeff, Sig])

# Define variable names for the results table
# Names for lagged returns
rString = [f"R_t-{lag}" for lag in range(1, maxlag + 1)]
# Names for weekday dummies
varnames = ["Mon", "Tue", "Wed", "Thu", "TaxDays"]
# Names for event dummies
eString = [f"Event +{fd_idx + 1}" for fd_idx in range(fD)]
# All variable names, including the constant (added by statsmodels)
var = ["Constant"] + rString + varnames + eString

# Create a pandas DataFrame to display the regression results
tab = pd.DataFrame(res, columns=var, index=['Coefficient', 't-statistic'])

# --- Display and Save Results ---
print("Regression Analysis Results:")
print(tab)

# Write the results table to an Excel file based on the America flag
if America == 0:
    tab.to_excel('RegressionAnalysis_All.xlsx')
else:
    tab.to_excel('RegressionAnalysis_America.xlsx')


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
