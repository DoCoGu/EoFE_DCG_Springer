# CAPM_3F.jl
# This script tests the Fama-French three-factor model for a set of assets
# by regressing excess asset returns on the market, SMB, and HML factors in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, GLM, Statistics

# Set working directory to current file location
cd(@__DIR__)

# === Data Loading and Preparation ===

# Load asset returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Load Fama-French Factors from an Excel file
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

# Extract Fama-French factors (Mkt-RF, SMB, HML) for the regression
# Select columns 2 to 4 (Mkt-RF, SMB, HML) from the second row (index 2) onwards
# Convert to a Matrix and divide by 100 to convert from percentage to decimal
Factors = Matrix(Factors[2:end, 2:4]) ./ 100  # size: T x 3 (Observations x Factors)

# === Model Estimation (Fama-French 3-Factor Model) ===

# Get the size of the Excess Returns matrix
nS, mS = size(ExRet)  # nS = number of observations, mS = number of assets

# Initialize arrays to store regression results for the 3-Factor model
# Use Array{Union{Missing, Float64}} to allow for potential missing values or NaNs
alpha = Array{Union{Missing, Float64}}(undef, 2, mS)  # Alpha coefficient and its t-statistic
beta_1 = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for Market factor and its t-statistic
beta_2 = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for SMB factor and its t-statistic
beta_3 = Array{Union{Missing, Float64}}(undef, 2, mS) # Beta for HML factor and its t-statistic
R2 = Array{Union{Missing, Float64}}(undef, mS)    # Adjusted R-squared

# Loop through each asset to perform the Fama-French 3-Factor regression
for i in 1:mS
    # Create a temporary DataFrame for the current asset's excess returns and the Fama-French factors
    # Use the columns from the Factors matrix directly
    df = DataFrame(Mkt=Factors[:, 1], ExRet=ExRet[:, i], SMB=Factors[:, 2], HML=Factors[:, 3])

    # Fit linear model using GLM: ExRet_i = alpha_i + beta1_i*Mkt + beta2_i*SMB + beta3_i*HML + epsilon_i
    model = lm(@formula(ExRet ~ Mkt + SMB + HML), df)

    # Extract regression results
    coefs = coef(model)         # Coefficients (alpha, beta_Mkt, beta_SMB, beta_HML)
    ctable = coeftable(model)   # Coefficient table containing statistics
    tstat = ctable.cols[3]      # Extract the t-statistics (assuming t-values are in the 4th column, index 3)

    # Store the regression results for the current asset
    alpha[1, i] = coefs[1]    # Alpha coefficient (intercept)
    alpha[2, i] = tstat[1]    # Alpha t-statistic
    beta_1[1, i] = coefs[2]   # Beta coefficient for Market
    beta_1[2, i] = tstat[2]   # Beta t-statistic for Market
    beta_2[1, i] = coefs[3]   # Beta coefficient for SMB
    beta_2[2, i] = tstat[3]   # Beta t-statistic for SMB
    beta_3[1, i] = coefs[4]   # Beta coefficient for HML
    beta_3[2, i] = tstat[4]   # Beta t-statistic for HML

    R2[i] = adjr2(model)      # Adjusted R-squared of the model
end

# === Format Results for Excel ===

# Prepare the results matrix for display and saving
# Vertically concatenate alpha, beta1, beta2, beta3, and R2 arrays
result_matrix = round.([alpha; beta_1; beta_2; beta_3; R2'], digits=5) # Transpose R2 for concatenation

# Define column names (asset names) from the original Returns DataFrame
col_names = names(Ret)[2:end]
# Define row names (statistics) for the results table
row_names = [
    "alpha", "(alpha t-stat)",
    "beta_mkt", "beta_mkt (t-stat)",
    "beta_smb", "(beta_smb t-stat)",
    "beta_hml", "(beta_hml t-stat)",
    "Adj. R2"
]

# Create a DataFrame from the results matrix
Results = DataFrame(result_matrix', :auto) # Transpose result_matrix and let DataFrame infer column names

# Rename the columns of the Results DataFrame to the defined row names
rename!(Results, Symbol.(row_names))

# Transpose the DataFrame to have statistics as rows and assets as columns
Results = permutedims(Results)

# Rename the columns to the defined column names (asset names)
rename!(Results, col_names)

# Add a new column 'Statistic' with the defined row names
Results = hcat(DataFrame(Statistic = row_names), Results)

# === Write Results to Excel ===
# Write the results table to an Excel file
XLSX.writetable("CAPM_3F.xlsx", Results; overwrite=true) # Specify overwrite=true to replace if file exists

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
