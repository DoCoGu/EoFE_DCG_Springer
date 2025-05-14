# CAPM.jl
# This script tests the single-factor Capital Asset Pricing Model (CAPM)
# for a set of assets by regressing excess asset returns on the market
# excess return in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using DataFrames, CSV, GLM, XLSX, Statistics

# Set working directory to current file location
cd(@__DIR__)

# === Data Loading and Preparation ===

# Import Stock Returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Import Fama-French Factors from an Excel file
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Extract Market Return and Risk-Free Rate from the Factors DataFrame
# Select rows from the second row (index 2) to the end
# Select the 2nd column (Market) and 5th column (Risk-Free Rate)
# Convert from percentage to decimal by dividing by 100
Mkt = Factors[2:end, 2] ./ 100  # Market return
Rf = Factors[2:end, 5] ./ 100  # Risk-free rate

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert the asset returns (columns 2 to end of Ret) to a Matrix
# Subtract the Risk-Free Rate (Rf) from each asset's return for each period
ExRet = Matrix(Ret[:, 2:end]) .- Rf

# === Testing the CAPM (One-Factor Model) ===

# Get the size of the Excess Returns matrix
nS, mS = size(ExRet)  # nS = number of observations, mS = number of assets

# Initialize arrays to store regression results
# Use Array{Union{Missing, Float64}} to allow for potential missing values or NaNs
alpha_S = Array{Union{Missing, Float64}}(undef, 2, mS)  # Alpha coefficient and its t-statistic
beta_S = Array{Union{Missing, Float64}}(undef, 2, mS)   # Beta coefficient and its t-statistic
R2_S = Array{Union{Missing, Float64}}(undef, mS)    # Adjusted R-squared

# Loop through each asset to perform the CAPM regression
for i in 1:mS
    # Create a temporary DataFrame for the current asset's excess returns and the market return
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i])

    # Fit the linear regression model using GLM
    # Model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
    model = lm(@formula(ExRet ~ Mkt), df)

    # Extract regression results
    coefs = coef(model)         # Coefficients (alpha and beta)
    ctable = coeftable(model)   # Coefficient table containing statistics
    tstat = ctable.cols[3]      # Extract the t-statistics (assuming t-values are in the 4th column, index 3)

    # Store the regression results for the current asset
    alpha_S[1, i] = coefs[1]    # Alpha coefficient (intercept)
    alpha_S[2, i] = tstat[1]    # Alpha t-statistic
    beta_S[1, i] = coefs[2]     # Beta coefficient (for Mkt)
    beta_S[2, i] = tstat[2]     # Beta t-statistic (for Mkt)
    R2_S[i] = adjr2(model)      # Adjusted R-squared of the model
end

# === Display and Save Results ===

# Prepare the results matrix for display and saving
# Vertically concatenate alpha, beta, and R2 arrays
result_matrix = round.([alpha_S; beta_S; R2_S'], digits=5) # Transpose R2_S for concatenation

# Define row and column names for the results table
row_names = ["alpha", "(alpha t-stat)", "beta", "(beta t-stat)", "Adj. R2"]
col_names = names(Ret)[2:end] # Use original asset names from the Returns DataFrame

# Create a DataFrame from the results matrix
Results = DataFrame(result_matrix', :auto) # Transpose result_matrix and let DataFrame infer column names

# Rename the columns of the Results DataFrame to the defined row names
rename!(Results, Symbol.(row_names))

# Transpose the DataFrame to have statistics as rows and assets as columns
Results = permutedims(Results)

# Rename the columns to the defined column names (asset names)
rename!(Results, col_names)

# Rename the columns again using Symbols (this line seems redundant if the previous rename worked)
# rename!(Results, Symbol.(col_names))

# Add a new column 'Statistic' with the defined row names
Results = hcat(DataFrame(Statistic = row_names), Results)

# Write the results table to an Excel file
XLSX.writetable("CAPM_Stock.xlsx", Results; overwrite=true) # Specify overwrite=true to replace if file exists

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
