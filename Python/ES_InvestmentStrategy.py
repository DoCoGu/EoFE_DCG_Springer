# ES_InvestmentStrategy.py
# This script analyzes different investment strategies, including buy-and-hold
# and an event-driven strategy based on aviation disasters, in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# --- Setting Parameters ---
cutoff = 150  # Cutoff value for filtering events (based on number of casualties)
America = 1  # Flag to filter events by Zone (1 for America, 0 for All)

# --- Data Loading and Preparation ---
# Import industry portfolio returns (assuming it includes Transport industry)
Ret_data = pd.read_excel("../Data/IndustryPortfolios.xlsx") # Renamed to avoid conflict with Ret variable later
# Import daily Fama-French Factors (specifically for Risk-Free Rate)
Rf_data = pd.read_excel("../Data/FF_FactorsDaily.xlsx")
# Extract Risk-Free Rate (RF) from the 'RF' column
Rf = Rf_data['RF'].values

# Extract dates from the returns data and convert to datetime objects
DatesReturn = pd.to_datetime(Ret_data.iloc[:, 0], format='%Y%m%d')

# Import Events data
Events_data = pd.read_excel('../Data/Events.xlsx') # Renamed to avoid conflict with Events variable later

# Filter events based on the America flag and casualty cutoff
if America == 0:
    # Filter for all events above the cutoff (assuming casualties are in column 12, event date in column 3)
    Events = Events_data[Events_data.iloc[:, 12] > cutoff].iloc[:, 3]
else:
    # Filter for events in Zone 1 (America, assuming Zone is in the last column) above the cutoff
    Events = Events_data[(Events_data.iloc[:, 12] > cutoff) & (Events_data["Zone"] == 1)].iloc[:, 3]

# Extract Transport industry returns (assuming it's in the 14th column, index 13)
Ret = Ret_data.iloc[:, 13].values

# Convert event dates to datetime objects (handling potential NaT from filtering)
# Use .values.flatten() to ensure a 1D numpy array before converting to datetime
Events = pd.to_datetime(Events.values.flatten(), format='%d-%m-%Y', errors='coerce')
# Remove any dates that failed to convert (resulted in NaT)
Events = Events.dropna()

# --- Create Event Dummy Variables ---
# Check if each date in DatesReturn is present in the filtered Events dates
Match = np.isin(DatesReturn, Events)

# Initialize event dummy matrix. Add extra rows initially to handle shifting for future events.
# Shape: (number of return dates + 3, 3) - Assuming 3 days after the event are considered
EDummy = np.zeros((len(DatesReturn) + 3, 3))

# Create dummies for event days and subsequent days (up to 3 days)
# Loop through each date index in DatesReturn
for i in range(len(DatesReturn)):
    # If the current date is an event date (found in the filtered Events)
    if Match[i]:
        # Set dummies for the next 3 days (indices 0, 1, 2 relative to the EDummy columns)
        # Check if the index (i + 1 + fd_idx) is within the bounds of the EDummy matrix
        if i + 1 < len(DatesReturn): # Check for day + 1
            EDummy[i + 1, 0] = 1
        if i + 2 < len(DatesReturn): # Check for day + 2
            EDummy[i + 2, 1] = 1
        if i + 3 < len(DatesReturn): # Check for day + 3
            EDummy[i + 3, 2] = 1

# Remove the first 3 rows of EDummy as they cannot have event dummies referring to previous events
# This aligns EDummy with the length of DatesReturn
EDummy = EDummy[3:, :]

# Create a dummy variable for tax days (Jan 1-5) - Note: This dummy is created but not used in the strategies below.
T = np.zeros(len(Ret))
dday = DatesReturn.dt.day
dmonth = DatesReturn.dt.month
T[(dmonth == 1) & (dday >= 1) & (dday <= 5)] = 1

# --- Investment Strategies Analysis ---

# Buy and Hold strategy: Cumulative excess return over time
# Calculate excess returns (Asset Return - Risk-Free Rate)
ExcessRet = Ret - Rf
# Calculate the cumulative sum of excess returns
BuyandHold = np.cumsum(ExcessRet)

# Aviation disasters event-driven strategy
# This strategy seems to modify the excess return on the day after an aviation disaster
# (where EDummy[:, 0] == 1) to be (Rf - Ret) instead of (Ret - Rf).
# This implies a specific action taken on those days.
AviationStrategy = ExcessRet.copy() # Start with excess returns
# Identify the indices where the dummy for the day after an event is 1 (EDummy[:, 0])
aviation_idx = (EDummy[:, 0] == 1)
# On these specific days, replace the excess return calculation as per the original code: Rf - Ret
# Ensure Rf and Ret are aligned for these indices
AviationStrategy[aviation_idx] = Rf[aviation_idx] - Ret[aviation_idx]

# Calculate the cumulative return for the aviation strategy
AviationStrategy = np.cumsum(AviationStrategy)

# --- Plotting the Strategies ---
# Create a figure with specified size (matching original code)
Inv_plot = plt.figure(1, figsize=(15, 10))

# Plot the Buy and Hold strategy cumulative returns
# Plot against the full DatesReturn
plt.plot(DatesReturn, BuyandHold, '-k', linewidth=1.5, label="Buy and Hold")
# Plot the Aviation Disaster Strategy cumulative returns
# Plot against the full DatesReturn (as in the original code, despite EDummy trimming)
# Note: This means the first 3 values of AviationStrategy correspond to the first 3 dates in DatesReturn,
# even though the EDummy values for these dates (before trimming) were zero.
plt.plot(DatesReturn, AviationStrategy, '--b', linewidth=1.5, label="Aviation disasters Strategy")

# Configure plot appearance
plt.legend(fontsize=16, loc="best") # Add legend with specified font size and location
plt.ylabel('Cumulative return', fontsize=18) # Set y-axis label with specified font size
plt.xlabel('Date', fontsize=18) # Set x-axis label with specified font size
plt.tight_layout() # Adjust layout to prevent labels overlapping

# Save the plot as an EPS file with specified format (matching original code)
plt.savefig('ES_InvStrategy.eps', format='eps')

# Show the plot (as in the original code)
plt.show()


# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
