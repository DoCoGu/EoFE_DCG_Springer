# ES_RegAnalysis.jl
# This script performs regression analysis on industry portfolio returns,
# including event studies and calendar effects (weekday and tax days) in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, Dates, GLM, StatsModels, ShiftedArrays, Statistics # Added Statistics for mean/std if needed, though not explicitly used in this version

# Set working directory to current file location
cd(@__DIR__)

# --- Settings ---
maxlag = 3       # Maximum lag for lagged returns in the regression
cutoff = 150     # Cutoff value for filtering events (based on number of casualties)
fD = 3           # Number of future days to create event dummies for
America = 0      # Flag to filter events by Zone (1 for America, 0 for All)

# --- Import Data ---

# Read industry portfolio returns and dates from an Excel file
# Read the "Portfolios" sheet from the Excel file into a DataFrame
Portfolios = DataFrame(XLSX.readtable("../Data/IndustryPortfolios.xlsx", "Portfolios"))
# Extract dates from the first column and convert to Date objects (assuming format "yyyymmdd")
DatesReturn = Date.(string.(Portfolios[:, 1]), dateformat"yyyymmdd")

# Read events data from an Excel file
# Read the "Events" sheet from the Excel file into a DataFrame
Events = DataFrame(XLSX.readtable("../Data/Events.xlsx", "Events"))

# --- Filter Events ---
# Filter events based on the America flag and casualty cutoff
# Assuming casualties are in column 13 and Zone is in the column named "Zone"
if America == 0
    # Filter for all events above the cutoff
    Events = Events[Events[:, 13] .> cutoff, 4] # Select the 4th column (event date)
else
    # Filter for events in Zone 1 (America) above the cutoff
    Events = Events[(Events[:, 13] .> cutoff) .& (Events[!, "Zone"] .== 1), 4] # Select the 4th column (event date)
end

# Convert filtered event dates to Date objects (assuming format "dd-mm-yyyy")
Events = Date.(Events[:, 1], dateformat"dd-mm-yyyy")

# Get the Transport industry return series (assuming it's in column 14)
Ret = Portfolios[:, 14]
nObs = size(Ret, 1) # Number of observations

# --- Weekday Dummies (Mon to Thu) ---
# Initialize array for weekday dummies (Monday=1, Tuesday=2, ..., Thursday=4)
# We create dummies for Monday (dayofweek 1) to Thursday (dayofweek 4)
dow = zeros(Int, nObs, 4)
for d = 1:4 # Loop for Monday (1) to Thursday (4)
    for i = 1:nObs # Loop through each observation date
        # Check if the weekday of the current date matches the current dummy day (d)
        if dayofweek(DatesReturn[i]) == d
            dow[i, d] = 1 # Set dummy to 1 if it's that weekday
        end
    end
end

# --- Event Dummies ---
# Initialize event dummy matrix. Adding extra rows initially to handle shifting.
# Shape: (number of return dates + fD, fD)
EDummy = zeros(Int, nObs + fD, fD)

# Find the row index in DatesReturn that corresponds to each event date
match_idx = findall(x -> x in Events, DatesReturn)

# Create dummies for event days and subsequent days (up to fD)
# The dummy EDummy[i+d, d] is 1 if day i is an event day and we are looking at day i+d.
for i in match_idx # Loop through indices of event days in DatesReturn
    for d = 1:fD # Loop through future days relative to the event day (1 to fD)
        # Ensure the index (i + d) is within the bounds of the extended EDummy array
        if i + d <= nObs + fD # Note: The original code used nObs here, but EDummy is larger
            EDummy[i + d, d] = 1
        end
    end
end
# Remove the first fD rows as they cannot have event dummies referring to previous events
EDummy = EDummy[fD+1:end, :] # Adjusting index for 1-based Julia arrays

# --- Tax Day Dummy (Jan 1–5) ---
# Initialize array for tax day dummy
T = zeros(Int, nObs)
# Extract day and month from return dates
dday = day.(DatesReturn)
dmonth = month.(DatesReturn)
# Set dummy to 1 if the date is in January (month 1) and between the 1st and 5th day
T[(dmonth .== 1) .& (dday .>= 1) .& (dday .<= 5)] .= 1 # Use broadcasting (.==, .>=, .<=, .=)

# --- Lagged Returns ---
# Create lagged returns matrix using ShiftedArrays.lag
# hcat combines the lagged columns horizontally
RetLag = hcat([lag(Ret, l) for l in 1:maxlag]...)
# Replace initial NaNs introduced by lag with zeros (or handle as missing data depending on GLM)
# GLM.jl can handle missing data, so replacing with NaN might be more appropriate if not trimming.
# If trimming is intended later, setting to 0 might be okay, but be mindful of its impact on regression.
# Let's keep it as per original code, assuming trimming handles it.
RetLag[1:maxlag, :] .= 0

# --- Design Matrix ---
# Combine all independent variables into a single matrix X
# Convert to Float64 Matrix for regression
X = convert(Matrix{Float64}, hcat(RetLag, dow, T, EDummy[:, 1:fD])) # Select columns 1 to fD of EDummy

# Define variable names for the regression predictors
# Names for lagged returns
rStrings = ["R_{t-$(i)}" for i in 1:maxlag]
# Names for weekday dummies
varnames_dow = ["Mon", "Tue", "Wed", "Thu"] # Corrected names for the 4 weekday dummies
# Name for tax day dummy
varnames_tax = ["TaxDays"]
# Names for event dummies
eStrings = ["Event +$(i)" for i in 1:fD]
# Combine all variable names
all_varnames = vcat(rStrings, varnames_dow, varnames_tax, eStrings)

# Create a DataFrame with named predictors and the dependent variable
df = DataFrame(X, Symbol.(all_varnames)) # Use all_varnames for predictor names
df.Returns = convert(Vector{Float64}, Ret) # Add the dependent variable (Returns)

# --- Regression ---
# Define the regression formula using StatsModels
# Returns ~ predictor1 + predictor2 + ... + predictorN
# Use `+ sum(Term.(Symbol.(all_varnames)))` to include all predictors dynamically
# Note: GLM's lm automatically includes an intercept unless `0 +` is used in the formula.
# The original MATLAB code used `fitlm(X, Ret(:,1))` which includes an intercept.
# Let's define the formula to include an intercept explicitly.
formula = Term(:Returns) ~ sum(Term.(Symbol.(all_varnames))) # Add `+ 1` for intercept

# Fit the linear regression model using GLM
# Select rows from index `maxlag + 1` to the end to remove rows with initial lagged NaNs/zeros
model = lm(formula, df[maxlag+1:end, :])

# --- Results ---
# Extract coefficients and standard errors
coeff = round.(coef(model), digits=5) # Coefficients
se = stderror(model) # Standard errors
# Calculate t-statistics (handle potential division by zero)
Sig = coeff ./ se # Element-wise division

# Combine coefficients and t-statistics into a single matrix
res = vcat(coeff', Sig') # Vertically concatenate transposed coefficient and t-statistic vectors

# --- Variable Names for Output Table ---
# Define variable names including the constant
output_varnames = ["Constant"; all_varnames]

# --- Output Table ---
# Create a DataFrame to display the regression results
# Reshape the results matrix (res) to have 2 rows (Coefficient, t-Stat) and columns for each variable
results_df = DataFrame(reshape(res, 2, :), Symbol.(output_varnames)) # Use output_varnames for column names
# Add a column for the row labels ('Coefficient', 't-Stat')
results_df[!, :Row] = ["Coefficient", "t-Stat"]
# Reorder columns to have 'Row' first
results_df = select(results_df, :Row, Not(:Row))

# --- Save Output ---
# Define the output file name based on the America flag
output_file = America == 0 ? "RegressionAnalysis_All.xlsx" : "RegressionAnalysis_America.xlsx"
# Write the results DataFrame to an Excel file, specifying the sheet name as "Results"
XLSX.writetable(output_file, "Results" => results_df; overwrite=true) # Specify overwrite=true

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
