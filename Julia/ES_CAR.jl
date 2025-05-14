# ES_CAR.jl
# This script performs an event study analysis to calculate Cumulative
# Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR)
# around specific events, using a CAPM model for expected returns in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, CSV, DataFrames, Dates, Statistics, LinearAlgebra, GLM, Plots # Added LinearAlgebra for matrix operations, Plots for plotting

# Set working directory to current file location
cd(@__DIR__)

# --- Settings ---
America = 0      # Flag: 1 to use only aviation disasters in America, 0 for all
cutoff = 150     # Cutoff for number of casualties to filter events

# Define the Cumulative Abnormal Return (CAR) window relative to event day (day 0)
start_CAR = -1   # Start of the window (e.g., -1 means 1 day before event)
end_CAR = 5      # End of the window (e.g., 5 means 5 days after event)

# Define the CAPM estimation window relative to event day (day 0)
start_CAPM = -250 # Start of the window
end_CAPM = -50    # End of the window

# --- Data Loading and Preparation ---
# Import industry portfolio returns (assuming it includes Transport industry)
# Read the "Portfolios" sheet from the Excel file into a DataFrame
Portfolios = DataFrame(XLSX.readtable("../Data/IndustryPortfolios.xlsx", "Portfolios"))

# Import daily Fama-French Factors
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_FactorsDaily.xlsx", "Factors"))

# Import Events data
# Read the "Events" sheet from the Excel file into a DataFrame
Events = DataFrame(XLSX.readtable("../Data/Events.xlsx", "Events"))

# --- Filter Events ---
# Filter events based on the America flag and casualty cutoff
# Assuming casualties are in column 13 and Zone is in the column named "Zone"
if America == 0
    # Filter for all events above the cutoff
    Events = Events[(Events[:, 13] .> cutoff), 4] # Select the 4th column (event date)
else
    # Filter for events in Zone 1 (America) above the cutoff
    Events = Events[(Events[:, 13] .> cutoff) .& (Events[:, "Zone"] .== 1), 4] # Select the 4th column (event date)
end

# --- Convert Returns Data ---
# Extract Market Excess Return (MktRet) and Risk-Free Rate (Rf) from Factors
# Select the 2nd column (Mkt-RF) and 5th column (RF) and divide by 100 to convert from percentage to decimal
MktRet = Factors[:, 2] ./ 100
Rf = Factors[:, 5] ./ 100

# Extract Transport industry returns (assuming it's in the 14th column)
# Select the 14th column and divide by 100 to convert from percentage to decimal
Ret = Portfolios[:, 14] ./ 100
# Calculate Excess Returns for the Transport industry (Asset Return - Risk-Free Rate)
ExRet = Ret .- Rf
# Reshape ExRet to a 2D matrix (n x 1) for consistency with matrix operations
ExRet = reshape(ExRet, :, 1)

# --- Convert Dates ---
# Convert event dates to Date objects (assuming format "dd-mm-yyyy")
Events_dates = Date.(Events[:, 1], dateformat"dd-mm-yyyy")
# Convert returns dates to Date objects (assuming format "yyyymmdd")
DatesReturn = Date.(string.(Factors[:, 1]), dateformat"yyyymmdd")

# Get the size of the Excess Returns matrix
n, nS = size(ExRet) # n = number of observations, nS = number of stocks (should be 1)

# --- Match Events to Return Dates ---
# Find the row index in DatesReturn that corresponds to each event date
Dates_vec = zeros(Int, length(Events_dates))
for i in 1:length(Events_dates)
    event_date = Events_dates[i]
    # Find the index of the event date in the returns dates
    idx_date = findfirst(x -> x == event_date, DatesReturn)
    # If the exact date is not found, search for the next available date
    while isnothing(idx_date)
        event_date += Day(1) # Move to the next day
        idx_date = findfirst(x -> x == event_date, DatesReturn)
        # Note: This loop assumes the event date or a subsequent date will eventually be found
        # within the DatesReturn range. Consider adding a safeguard against infinite loops
        # if an event date is after the last return date.
    end
    Dates_vec[i] = idx_date # Store the found index
end

# Filter out event dates where the estimation window would go before the start of the data
# Keep only event dates where the start of the CAPM estimation window (Dates_vec[i] + start_CAPM) is >= 1
Dates_event = Dates_vec[Dates_vec .> abs(start_CAPM)] # Assuming start_CAPM is negative
Dates_event = unique(Dates_event) # Ensure unique event dates
nE = length(Dates_event) # Number of valid events

# --- Event Study: Estimate CAPM Alpha and Beta ---
# Initialize arrays to store CAPM parameters (alpha and beta) for each valid event
alpha = fill(NaN, nE, nS)
beta = fill(NaN, nE, nS)

# Loop through each valid event to estimate CAPM parameters in the estimation window
for i in 1:nE
    # Define the indices for the estimation window relative to the event date index
    estimation_window_start_idx = Dates_event[i] + start_CAPM
    estimation_window_end_idx = Dates_event[i] + end_CAPM

    # Ensure the estimation window is within the bounds of the data (1 to n)
    if estimation_window_start_idx >= 1 && estimation_window_end_idx <= n
        # Extract the data for the current estimation window
        X = MktRet[estimation_window_start_idx : estimation_window_end_idx] # Market Excess Return
        Y = ExRet[estimation_window_start_idx : estimation_window_end_idx, :] # Asset Excess Return (for current asset k)

        # Create a temporary DataFrame for the regression
        df_n = DataFrame(X = X, Y = Y[:, 1]) # Assuming only one asset (nS=1)

        # Fit CAPM model using GLM: Y = alpha + beta * X + epsilon
        model = lm(@formula(Y ~ X), df_n)

        # Store the estimated alpha and beta for the current event and asset
        # coef(model)[1] is the intercept (alpha), coef(model)[2] is the coefficient for X (beta)
        alpha[i, 1] = coef(model)[1] # Assuming nS=1
        beta[i, 1] = coef(model)[2]  # Assuming nS=1
    end
end

# --- Calculate Predicted Returns ---
# Define the length of the CAR window
len_CAR = abs(start_CAR) + abs(end_CAR) + 1

# Initialize array to store predicted returns around events
# Shape: (length of CAR window, number of events, number of stocks)
PredRet = fill(NaN, len_CAR, nE, nS)

# Calculate predicted returns for each valid event and each day in the CAR window
for t in 1:len_CAR # Loop through days in the CAR window (relative index)
    for i in 1:nE # Loop through valid events
        for k in 1:nS # Loop through stocks (should be 1)
            # Calculate the index in the returns data for the current day in the CAR window
            # Dates_event[i] is the index of the event day (day 0)
            # start_CAR is the offset to the start of the CAR window
            # t is the current position within the CAR window (1 to len_CAR)
            prediction_date_idx = Dates_event[i] + start_CAR + (t - 1) # Adjusting index for 1-based Julia arrays

            # Ensure the date for prediction is within the bounds of the data (1 to n)
            if prediction_date_idx >= 1 && prediction_date_idx <= n
                 # Check if CAPM parameters were successfully estimated for this event
                if !isnan(alpha[i, k]) && !isnan(beta[i, k])
                    # Predict return using estimated CAPM parameters and market return on that day
                    PredRet[t, i, k] = alpha[i, k] + beta[i, k] * MktRet[prediction_date_idx]
                end
            end
        end
    end
end

# --- Get Observed Returns ---
# Preallocate array for observed returns around events
# Shape: (length of CAR window, number of events, number of stocks)
ObsRet_agg = fill(NaN, len_CAR, nE, nS)

# Get observed returns for each valid event and each day in the CAR window
for i in 1:nE # Loop through valid events
    for t in 1:len_CAR # Loop through days in the CAR window (relative index)
        for k in 1:nS # Loop through stocks (should be 1)
            # Calculate the index in the returns data for the current day in the CAR window
            observation_date_idx = Dates_event[i] + start_CAR + (t - 1) # Adjusting index for 1-based Julia arrays

            # Ensure the date for observation is within the bounds of the data (1 to n)
            if observation_date_idx >= 1 && observation_date_idx <= n
                ObsRet_agg[t, i, k] = ExRet[observation_date_idx, k] # Get observed excess return
            end
        end
    end
end

# --- Calculate Abnormal and Cumulative Returns ---
# Get abnormal returns (Observed return - Predicted return)
AbnRet = ObsRet_agg .- PredRet

# Get cumulative abnormal returns (sum of abnormal returns over the CAR window)
# mapslices applies a function along a specified dimension.
# cumsum(skipmissing(x)) calculates cumulative sum ignoring missing values within each slice.
# dims=1 applies the function along the first dimension (time).
CAR = mapslices(x -> cumsum(collect(skipmissing(x))), AbnRet, dims=1) # Use collect to handle SkipMissing iterator

# Get cumulative average abnormal returns (average of CARs across all events)
# mean(..., dims=2) calculates the mean along the second dimension (events).
# [:] flattens the resulting 1x1x1 array (if nS=1) into a vector.
CAAR = mean(CAR, dims=2)[:]
CAAR .*= 100 # Convert to percentage

# --- Plotting CAAR ---
# Get dates vector for the plot (relative to event day 0)
date = collect(start_CAR:end_CAR) # Create a range of integers for the x-axis

# Create a zero vector for the horizontal line at zero
zero_vec = zeros(length(date))

# Create the plot for Cumulative Average Abnormal Returns (CAAR)
p = plot(date, CAAR, lw=2, label="CAAR", color=:blue) # Added label and color

# Add a horizontal line at zero for reference
plot!(p, date, zero_vec, lw=2, linestyle=:dash, color=:black, label="") # Plot onto plot p

# Add a vertical line at day 0 (the event day) for reference
vline!(p, [0], lw=2, color=:black, linestyle=:solid, label="") # Plot onto plot p

# Configure plot appearance
xlabel!("Days Relative to Events") # X-axis label
ylabel!("CAAR (%)") # Y-axis label (in percentage)
xlims!(start_CAR, end_CAR) # Set x-axis limits
# Set y-axis limits with a small buffer
ylims!(minimum(CAAR) - abs(minimum(CAAR)*0.1), maximum(CAAR) + abs(maximum(CAAR)*0.1)) # Adjusted buffer calculation

title!("Cumulative Average Abnormal Returns (CAAR)") # Plot title
plot!(p, legend=:topleft) # Legend position

# Save the plot as a PNG file, with filename based on the America setting
savefig(p, America == 0 ? "CAAR_All.png" : "CAAR_America.png") # Specify plot object and filename

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
