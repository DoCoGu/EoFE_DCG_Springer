# CS_EPU.jl
# This script performs Fama-MacBeth cross-sectional regression analysis,
# including Economic Policy Uncertainty (EPU) and consumption growth as factors,
# and calculates Newey-West and Shanken standard errors in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, Statistics, LinearAlgebra, GLM, CSV, CovarianceMatrices

# Set working directory to current file location
cd(@__DIR__)

# --- Helper Function: Newey-West Standard Errors ---
"""
    newey_west_se(X::Matrix, y::Vector, bandwidth::Int)

Calculates Newey-West (HAC) standard errors for OLS regression coefficients.

Parameters:
    X : Matrix
        Design matrix (independent variables).
    y : Vector
        Dependent variable.
    bandwidth : Int
        The number of lags to include in the HAC covariance estimation.

Returns:
    b : Vector
        OLS regression coefficients.
    se_nw : Vector
        Newey-West standard errors for the coefficients.
"""
function newey_west_se(X::Matrix, y::Vector, bandwidth::Int)
    n, k = size(X) # n = number of observations, k = number of coefficients
    b = X \ y # Calculate OLS coefficients using backslash operator (least squares)
    u = y .- X * b # Calculate residuals

    # Initialize the HAC covariance matrix (S)
    S = zeros(k, k)

    # Loop through lags from 0 to the specified bandwidth
    for l in 0:bandwidth
        # Calculate the Newey-West weight for the current lag
        w_l = 1.0
        if l > 0
            w_l -= l / (bandwidth + 1)
        end

        # Calculate the autocovariance matrix for the current lag (Gamma_l)
        Γ_l = zeros(k, k)
        # Iterate through observations to compute the sum for Gamma_l
        for t in (l + 1):n
            xt = X[t, :]' # Transpose the row vector to get a column vector
            xt_l = X[t - l, :] # Row vector from the lagged observation
            # Accumulate the product of lagged design matrices and residuals
            Γ_l = Γ_l .+ xt * xt_l * u[t] * u[t - l]
        end

        # Add the weighted autocovariance to the HAC covariance matrix S
        # For lag 0, add Γ_0. For lags > 0, add w_l * (Γ_l + Γ_l')
        S = S .+ w_l * (Γ_l + Γ_l')
        if l == 0
            # At lag 0, Γ_0 is added twice in the loop (as Γ_l and Γ_l'). Subtract one Γ_0.
            S = S .- Γ_l
        end
    end

    # Calculate the inverse of (X' * X)
    XtX_inv = inv(X' * X)
    # Calculate the variance-covariance matrix of the coefficients (V) using the sandwich estimator
    V = XtX_inv * S * XtX_inv
    # Return the OLS coefficients and the square root of the diagonal elements of V (standard errors)
    return b, sqrt.(diag(V))
end

# --- Data Import ---
# Load Portfolio Returns from an Excel file
Ret = DataFrame(XLSX.readtable("../Data/Portfolios.xlsx", "Portfolios"))
# Load Fama-French Factors from an Excel file
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Load Uncertainty data (including EPU) from an Excel file
Uncertainty = DataFrame(XLSX.readtable("../Data/Uncertainty.xlsx", "Uncertainty"))
# Extract EPU data and convert to a Float64 Vector
Uncertainty = convert(Vector{Float64}, Uncertainty[:, 3])

# Load Consumption data from an Excel file
Consumption = DataFrame(XLSX.readtable("../Data/Consumption.xlsx", "Consumption"))
# Extract Consumption data and convert to a Float64 Vector
Consumption = convert(Vector{Float64}, Consumption[:, 2])

# --- Data Preparation ---
# Compute EPU and Consumption Growth (percentage changes)
# EPU change: (Current EPU / Previous EPU) - 1
EPU = Uncertainty[2:end] ./ Uncertainty[1:end-1] .- 1
# Consumption change: (Current Consumption / Previous Consumption) - 1
C = Consumption[2:end] ./ Consumption[1:end-1] .- 1

# Extract Market Excess Return (Mkt) and Risk-Free Rate (Rf)
# Select the Market column (index 2) from the second row (index 2) onwards, convert to Float64 Vector, and divide by 100
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
# Select the Risk-Free Rate column (index 5) from all rows, convert to Float64 Vector, and divide by 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Calculate Excess Returns for portfolios (Portfolio Return - Risk-Free Rate)
# Convert the portfolio returns (columns 2 to end of Ret) from the second row (index 2) onwards to a Matrix and divide by 100
# Subtract the Risk-Free Rate (Rf) from the second row (index 2) onwards to match the length of other series
ExRet = (Matrix(Ret[2:end, 2:end]) / 100 .- Rf[2:end])

# Combine the factors (Market, Consumption Change, EPU Change) into a Matrix
Factors = hcat(Mkt, C, EPU)

# --- First Stage Regression (Time-Series Regressions) ---
# Get the size of the Excess Returns and Factors matrices
n1, n2 = size(ExRet) # n1 = number of observations, n2 = number of portfolios
nF = size(Factors, 2) # nF = number of factors

# Initialize arrays to store results from the first stage
# CoefAll stores the factor betas for each portfolio (excluding intercept)
CoefAll = fill(NaN, nF, n2)
# Res stores the residuals from each regression
Res = fill(NaN, n1, n2)

# Loop through each portfolio to perform time-series regression
for i in 1:n2
    # Create a temporary DataFrame for the current portfolio's excess returns and the factors
    df_n = DataFrame(Mkt=Factors[:, 1], ExRet=ExRet[:, i], C=Factors[:, 2], EPU=Factors[:, 3])

    # Fit linear model using GLM: ExRet_i = alpha_i + beta_Mkt*Mkt + beta_dC*dC + beta_EPU*EPU + epsilon_i
    model = lm(@formula(ExRet ~ Mkt + C + EPU), df_n)

    # Extract coefficients (alpha and betas)
    Coef = coef(model)
    # Store the factor betas (excluding the intercept, which is at index 1)
    CoefAll[:, i] = Coef[2:end]
    # Store the residuals
    Res[:, i] = residuals(model)
end

# Calculate the variance-covariance matrix of the residuals
VarCovErr = cov(Res)

# --- Second Stage Regression (Cross-Sectional Regressions) ---
# Calculate the average excess return for each portfolio
MeanRet = convert(Matrix{Float64}, mean(ExRet, dims=1)') # Calculate mean along dimension 1 and transpose to get a column vector

# Transpose the factor betas matrix from the first stage
Betas = convert(Matrix{Float64}, CoefAll') # Shape: (number of portfolios, number of factors)

# Create a DataFrame for the cross-sectional regression
df = DataFrame(hcat(Betas, MeanRet), :auto) # Combine Betas and MeanRet
rename!(df, ["Mkt", "C", "EPU", "MeanRet"]) # Rename columns

# Perform cross-sectional regression of average returns on factor betas
# Model: MeanRet_i = lambda_Mkt*beta_Mkt_i + lambda_dC*beta_dC_i + lambda_EPU*beta_EPU_i + eta_i
# Using `0 +` in the formula suppresses the intercept (assuming lambda_0 = 0)
model2 = lm(@formula(MeanRet ~ 0 + Mkt + C + EPU), df)

# Extract factor risk premia (Lambdas) and their standard errors
Lambda = coef(model2) # Factor risk premia
SE = stderror(model2) # Standard errors of Lambdas
Tstat = Lambda ./ SE # T-statistics (standard)

# --- Shanken Correction ---
# Calculate the covariance matrix of the factors
Sigma_f = cov(Factors)
# Calculate the inverse of (Betas' * Betas)
BtB_inv = inv(Betas' * Betas)
# Calculate the correction term for Shanken standard errors
correction = 1 + Lambda' * inv(Sigma_f) * Lambda
# Calculate the Shanken-corrected variance-covariance matrix of Lambdas
VarLam = (BtB_inv * Betas' * VarCovErr * Betas * BtB_inv * correction + Sigma_f) / n1
# Extract Shanken-corrected standard errors (square root of the diagonal)
SE_Shanken = sqrt.(diag(VarLam))
# Calculate Shanken-corrected T-statistics
Tstat_Shanken = Lambda ./ SE_Shanken

# --- Time-series of cross-section regressions (for Newey-West standard errors) ---
# Initialize matrix to store Lambdas from each cross-sectional regression
LambdaFull = fill(NaN, n1, nF) # Shape: (number of observations, number of factors)

# Loop through each time period (observation)
for j in 1:n1
    # Create a DataFrame for the cross-sectional regression for the current time period
    # Combine excess returns for the current period with the betas from the first stage
    df_t = DataFrame(hcat(ExRet[j, :], Betas), :auto)
    rename!(df_t, ["ExRet", "Mkt", "C", "EPU"]) # Rename columns

    # Fit cross-sectional model for the current time period (using betas from the first stage)
    # Model: ExRet_j_i = lambda_Mkt_j*beta_Mkt_i + lambda_dC_j*beta_dC_i + lambda_EPU_j*beta_EPU_i + eta_j_i
    # Using `0 +` in the formula suppresses the intercept
    model_j = lm(@formula(ExRet ~ 0 + Mkt + C + EPU), df_t)
    # Store the estimated Lambdas for this time period (transpose to get a row vector)
    LambdaFull[j, :] = coef(model_j)'
end

# Calculate the mean of the estimated Lambdas across all time periods
LambdaMean = mean(LambdaFull, dims=1)[:] # Calculate mean along dimension 1 and flatten to a vector

# Calculate HAC-corrected (Newey-West) standard errors
# This requires the time series of estimated Lambdas (LambdaFull)
nF = size(LambdaFull, 2) # Number of factors (same as number of Lambdas)
SE_NW = zeros(nF) # Initialize array for Newey-West standard errors

# Loop through each factor to calculate its HAC standard error using the custom function
for k in 1:nF
    y = LambdaFull[:, k] # Time series of estimated Lambdas for factor k
    X = ones(n1, 1) # Design matrix for the dummy regression (intercept only)
    bandwidth = 1 # Bandwidth for Newey-West (adjust as needed)

    # Compute HAC robust variance using the newey_west_se function
    # The function returns coefficients and standard errors
    b, SE_t = newey_west_se(X, y, bandwidth)

    # Extract the Newey-West standard error for the intercept (which is the mean of LambdaFull[:, k])
    SE_NW[k] = SE_t[1]
end

# Calculate Newey-West T-statistics
Tstat_NW = LambdaMean ./ SE_NW

# --- Results Table Formatting ---
# Get the names of the portfolios from the original Returns DataFrame
NamePort = names(Ret)[2:end]

# Create a DataFrame for the First Stage results (Factor Betas)
FirstStageReg = DataFrame(Betas, :auto) # Create DataFrame from Betas matrix
rename!(FirstStageReg, [:Mkt, :dC, :EPU]) # Rename columns to factor names
# Add a column for portfolio names
FirstStageReg.rowid = NamePort

# Create a DataFrame for the Second Stage results (Factor Risk Premia and T-stats)
SecondStage = DataFrame(
    Lambda = Lambda,
    tstat = Tstat,
    t_stat_HAC = Tstat_NW,
    t_stat_Shanken = Tstat_Shanken,
    rowindex = ["Mkt", "dC", "EPU"] # Row names for the statistics (factor names)
)
# Stack the DataFrame to long format
SecondStage = stack(SecondStage, Not(:rowindex))
# Unstack the DataFrame to wide format, with statistics as columns and factors as rows
SecondStage = unstack(SecondStage, :variable, :value)

# --- Write Results to XLSX ---
# Write the First Stage results table to a XLSX file
XLSX.writetable("FirstStage_EPU.xlsx", FirstStageReg; overwrite=true)
# Write the Second Stage results table to a XLSX file
XLSX.writetable("SecondStage_EPU.xlsx", SecondStage; overwrite=true)

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
