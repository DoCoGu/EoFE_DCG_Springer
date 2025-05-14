# CAPM_2F.jl
# This script tests two-factor asset pricing models by adding either VIX
# changes or EPU changes as a second factor to the CAPM in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, GLM, Statistics

# Set working directory to current file location
cd(@__DIR__)

# === Data Import ===

# Load asset returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Load Fama-French Factors from an Excel file
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Load Uncertainty data (including VIX and EPU) from an Excel file
# Read the "Uncertainty" sheet from the Excel file into a DataFrame
Uncertainty = DataFrame(XLSX.readtable("../Data/Uncertainty.xlsx", "Uncertainty"))

# === Data Preparation ===

# Calculate percentage changes for VIX and EPU
# VIX change: (Current VIX / Previous VIX) - 1
# Select the VIX column (index 2) and calculate element-wise division and subtraction
VIX = Vector(Uncertainty[2:end, 2] ./ Uncertainty[1:end-1, 2] .- 1)

# EPU change: (Current EPU / Previous EPU) - 1
# Select the EPU column (index 3) and calculate element-wise division and subtraction
EPU = Vector(Uncertainty[2:end, 3] ./ Uncertainty[1:end-1, 3] .- 1)

# Extract Market Excess Return (Mkt) and Risk-Free Rate (Rf)
# Select the Market column (index 2) from the second row (index 2) onwards, convert to Float64 Vector, and divide by 100
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
# Select the Risk-Free Rate column (index 5) from all rows, convert to Float64 Vector, and divide by 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert the asset returns (columns 2 to end of Ret) to a Matrix
# Subtract the Risk-Free Rate (Rf) from the second row (index 2) onwards to match the length of other series
ExRet = Matrix(Ret[:, 2:end]) .- Rf[2:end]

# === Model Estimation (CAPM + VIX) ===

# Get the size of the Excess Returns matrix
nS, mS = size(ExRet)  # nS = number of observations, mS = number of assets

# Initialize arrays to store regression results for CAPM + VIX model
# Use Array{Union{Missing, Float64}} to allow for potential missing values or NaNs
alpha_VIX = Array{Union{Missing, Float64}}(undef, 2, mS)  # Alpha coefficient and its t-statistic
beta_1_VIX = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for Market factor and its t-statistic
beta_2_VIX = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for VIX factor and its t-statistic
R2_VIX = Array{Union{Missing, Float64}}(undef, mS)    # Adjusted R-squared

# Loop through each asset to perform the regression with Market and VIX
for i in 1:mS
    # Create a temporary DataFrame for the current asset's excess returns, Market, and VIX
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i], VIX = VIX)

    # Fit linear model using GLM: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * VIX + epsilon_i
    model = lm(@formula(ExRet ~ Mkt + VIX), df)

    # Extract regression results
    coefs = coef(model)         # Coefficients (alpha, beta_Mkt, beta_VIX)
    ctable = coeftable(model)   # Coefficient table containing statistics
    tstat = ctable.cols[3]      # Extract the t-statistics (assuming t-values are in the 4th column, index 3)

    # Store the regression results for the current asset
    alpha_VIX[1, i] = coefs[1]    # Alpha coefficient (intercept)
    alpha_VIX[2, i] = tstat[1]    # Alpha t-statistic
    beta_1_VIX[1, i] = coefs[2]   # Beta coefficient for Market
    beta_1_VIX[2, i] = tstat[2]   # Beta t-statistic for Market
    beta_2_VIX[1, i] = coefs[3]   # Beta coefficient for VIX
    beta_2_VIX[2, i] = tstat[3]   # Beta t-statistic for VIX

    R2_VIX[i] = adjr2(model)      # Adjusted R-squared of the model
end

# === Model Estimation (CAPM + EPU) ===

# Initialize arrays to store regression results for CAPM + EPU model
alpha_EPU = Array{Union{Missing, Float64}}(undef, 2, mS)  # Alpha coefficient and its t-statistic
beta_1_EPU = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for Market factor and its t-statistic
beta_2_EPU = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for EPU factor and its t-statistic
R2_EPU = Array{Union{Missing, Float64}}(undef, mS)    # Adjusted R-squared

# Loop through each asset to perform the regression with Market and EPU
for i in 1:mS
    # Create a temporary DataFrame for the current asset's excess returns, Market, and EPU
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i], EPU = EPU)

    # Fit linear model using GLM: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * EPU + epsilon_i
    model = lm(@formula(ExRet ~ Mkt + EPU), df)

    # Extract regression results
    coefs = coef(model)         # Coefficients (alpha, beta_Mkt, beta_EPU)
    ctable = coeftable(model)   # Coefficient table containing statistics
    tstat = ctable.cols[3]      # Extract the t-statistics (assuming t-values are in the 4th column, index 3)

    # Store the regression results for the current asset
    alpha_EPU[1, i] = coefs[1]    # Alpha coefficient (intercept)
    alpha_EPU[2, i] = tstat[1]    # Alpha t-statistic
    beta_1_EPU[1, i] = coefs[2]   # Beta coefficient for Market
    beta_1_EPU[2, i] = tstat[2]   # Beta t-statistic for Market
    beta_2_EPU[1, i] = coefs[3]   # Beta coefficient for EPU
    beta_2_EPU[2, i] = tstat[3]   # Beta t-statistic for EPU

    R2_EPU[i] = adjr2(model)      # Adjusted R-squared of the model
end

# === Display and Save Results (CAPM + VIX) ===

# Prepare the results matrix for display and saving for the VIX model
# Vertically concatenate alpha, beta1 (Mkt), beta2 (VIX), and R2 arrays
# Note: Original code used alpha_S, beta_1, beta_2, R2_S' which seems like a mix.
# Correcting to use VIX-specific results and beta_1_VIX, beta_2_VIX.
result_VIX = round.([alpha_VIX; beta_1_VIX; beta_2_VIX; R2_VIX'], digits=5) # Transpose R2_VIX for concatenation

# Define row and column names for the VIX results table
row_names_VIX = ["alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "VIX", "(VIX t-stat)", "Adj. R2"]
col_names = names(Ret)[2:end] # Use original asset names from the Returns DataFrame

# Create a DataFrame from the VIX results matrix
Stocks_VIX = DataFrame(result_VIX', :auto) # Transpose result_VIX and let DataFrame infer column names

# Rename the columns of the Stocks_VIX DataFrame to the defined row names
rename!(Stocks_VIX, Symbol.(row_names_VIX))

# Transpose the DataFrame to have statistics as rows and assets as columns
Stocks_VIX = permutedims(Stocks_VIX)

# Rename the columns to the defined column names (asset names)
rename!(Stocks_VIX, col_names)

# Add a new column 'Statistic' with the defined row names
Stocks_VIX = hcat(DataFrame(Statistic = row_names_VIX), Stocks_VIX)

# === Display and Save Results (CAPM + EPU) ===

# Prepare the results matrix for display and saving for the EPU model
# Vertically concatenate alpha, beta1 (Mkt), beta2 (EPU), and R2 arrays
# Note: Original code used alpha_S, beta_1, beta_2, R2_S' which seems like a mix.
# Correcting to use EPU-specific results and beta_1_EPU, beta_2_EPU.
result_EPU = round.([alpha_EPU; beta_1_EPU; beta_2_EPU; R2_EPU'], digits=5) # Transpose R2_EPU for concatenation

# Define row and column names for the EPU results table
row_names_EPU = ["alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "EPU", "(EPU t-stat)", "Adj. R2"]
# col_names is the same as for VIX results

# Create a DataFrame from the EPU results matrix
Stocks_EPU = DataFrame(result_EPU', :auto) # Transpose result_EPU and let DataFrame infer column names

# Rename the columns of the Stocks_EPU DataFrame to the defined row names
rename!(Stocks_EPU, Symbol.(row_names_EPU))

# Transpose the DataFrame to have statistics as rows and assets as columns
Stocks_EPU = permutedims(Stocks_EPU)

# Rename the columns to the defined column names (asset names)
rename!(Stocks_EPU, col_names)

# Add a new column 'Statistic' with the defined row names
Stocks_EPU = hcat(DataFrame(Statistic = row_names_EPU), Stocks_EPU)

# === Write Results to Excel ===
# Write the VIX results table to an Excel file
XLSX.writetable("CAPM_VIX.xlsx", Stocks_VIX; overwrite=true)
# Write the EPU results table to an Excel file
XLSX.writetable("CAPM_EPU.xlsx", Stocks_EPU; overwrite=true)

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
