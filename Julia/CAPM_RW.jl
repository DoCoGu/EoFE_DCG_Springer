# CAPM_RW.jl
# This script implements a rolling-window regression to estimate time-varying
# CAPM alpha and beta for a set of assets and plots the results in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, Statistics, GLM, Plots, Dates

# Set working directory to current file location
cd(@__DIR__)

# --- Data Loading and Preparation ---

# Import Stock Returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Import Fama-French Factors from an Excel file
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Extract Market Excess Return (Mkt) and Risk-Free Rate (Rf)
# Select the Market column (index 2) from the second row (index 2) onwards, convert to Float64 Vector, and divide by 100
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
# Select the Risk-Free Rate column (index 5) from all rows, convert to Float64 Vector, and divide by 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert the asset returns (columns 2 to end of Ret) to a Matrix
# Subtract the Risk-Free Rate (Rf) from the second row (index 2) onwards to match the length of other series
ExRet = Matrix(Ret[:, 2:end]) .- Rf[2:end]

# --- Rolling Regression Setup ---

# Define the size of the rolling window
rw = 40 # Number of observations in each window

# Get the size of the Excess Returns matrix
nS, mS = size(ExRet)  # nS = total number of observations, mS = number of assets

# Initialize arrays to store rolling regression results
# Use fill(NaN, ...) to create arrays filled with NaN of the appropriate size
alpha = fill(NaN, nS - rw + 1, mS) # Alpha coefficient for each window and asset
beta = fill(NaN, nS - rw + 1, mS)  # Beta coefficient for each window and asset

# --- Perform Rolling Window Regression ---
# Outer loop iterates through the starting point of each window
for t in 1:(nS - rw + 1)
    # Inner loop iterates through each asset
    for i in 1:mS
        # Create a temporary DataFrame for the current window's data
        # Select the relevant rows for the current window (t to t+rw-1)
        df = DataFrame(ExRet=ExRet[t:t+rw-1, i], Mkt=Mkt[t:t+rw-1])

        # Fit the linear regression model for the current window and asset using GLM
        # Model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
        model = lm(@formula(ExRet ~ Mkt), df)

        # Extract and store the regression results (coefficients) for the current window
        alpha[t, i] = coef(model)[1] # Alpha coefficient (intercept)
        beta[t, i] = coef(model)[2]  # Beta coefficient (for Mkt)
    end
end

# --- Plotting Rolling Estimates ---

# Generate dates for the x-axis of the plots
# Select dates corresponding to the end of each rolling window (from row rw to the end of Ret)
dates = Ret[rw:end, 1]  # Assuming dates are in the first column of Ret
# Convert dates to Date objects
dates = Date.(dates)

# Generate a range of points for x-axis ticks based on the date range
# This creates 20 points evenly spaced across the duration of the rolling window results
start_date = dates[1]
end_date = dates[end]
day_range = LinRange(0, Dates.value(end_date - start_date), 20)
# Create Date objects for the plot ticks
Dates_for_plot = [start_date + Day(round(Int, d)) for d in day_range]

# Format the date labels for the x-axis ticks
date_labels = Dates.format.(DateTime.(Dates_for_plot), dateformat"u-yyyy") # e.g., "Jan-2000"

# Alpha plot
# Create the initial plot for the first asset's rolling alpha
plot(1:nS - rw + 1, alpha[:, 1], label=names(Ret)[2], lw=1.5, color=:red)
# Add rolling alpha for subsequent assets to the same plot
for i in 2:mS
    plot!(1:nS - rw + 1, alpha[:, i], label=names(Ret)[i+1], lw=1.5) # Use i+1 for correct column name index
end
# Add a horizontal line at zero for reference
plot!(1:nS - rw + 1, zeros(nS - rw + 1), ls=:dash, color=:black, label="")

# Configure x-axis ticks and labels
xticks!(range(1, stop=nS - rw + 1, length=20), date_labels)
plot!(xtickfont=font(8)) # Adjust font size for ticks
plot!(xrotation=45) # Rotate x-axis labels

# Set title and save the plot
title!("Rolling Alpha (CAPM)")
savefig("alpha_S.png") # Save as PNG file

# Beta plot
# Create the initial plot for the first asset's rolling beta
plot(1:nS - rw + 1, beta[:, 1], label=names(Ret)[2], lw=1.5, color=:red)
# Add rolling beta for subsequent assets to the same plot
for i in 2:mS
    plot!(1:nS - rw + 1, beta[:, i], label=names(Ret)[i+1], lw=1.5) # Use i+1 for correct column name index
end
# Add a horizontal line at one for reference (CAPM beta is expected to be 1 for the market)
plot!(1:nS - rw + 1, ones(nS - rw + 1), ls=:dash, color=:black, label="")

# Configure x-axis ticks and labels (same as alpha plot)
xticks!(range(1, stop=nS - rw + 1, length=20), date_labels)
plot!(xtickfont=font(8)) # Adjust font size for ticks
plot!(xrotation=45) # Rotate x-axis labels

# Set title and save the plot
title!("Rolling Beta (CAPM)")
savefig("Beta_S.png") # Save as PNG file


# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
