# BlackLitterman.jl
# This script implements the Black-Litterman model for portfolio allocation,
# blending market equilibrium returns with investor views, and plots the
# resulting efficient frontiers in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, GLM, Statistics, LinearAlgebra, Plots # Added LinearAlgebra for matrix operations, Plots for plotting

# Set working directory to current file location
cd(@__DIR__)

# Include the helper function for efficient frontier calculations
include("getuncef.jl")

# --- Data Loading and Preparation ---

# Import asset returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Import Fama-French Factors from an Excel file
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Add Market Excess Return and Risk-Free Rate columns to the returns DataFrame
# Select "Mkt-RF" and "RF" columns from the second row (index 2) onwards of the Factors DataFrame
# Convert from percentage to decimal by dividing by 100
Ret[!, "Mkt-RF"] = Factors[2:end, "Mkt-RF"] ./ 100
Ret[!, "rf"] = Factors[2:end, "RF"] ./ 100

# Calculate excess returns for assets (asset return - risk-free rate)
# Get the names of the asset columns (excluding Date, Mkt-RF, rf)
asset_columns = names(Ret)[2:end-2]
# Convert asset return columns to a Matrix and subtract the risk-free rate
exret = Matrix(Ret[:, asset_columns]) .- Ret.rf
# Extract Market Excess Return column
exmkt = Ret[!, "Mkt-RF"]

# --- Compute the Implied Equilibrium Returns (Pi) using CAPM ---
# This step estimates CAPM beta for each asset to derive equilibrium returns

T, n = size(exret) # T = number of observations, n = number of assets
beta = zeros(2, n) # Initialize matrix to store CAPM betas (intercept and slope)

# Loop through each asset to perform the CAPM regression (Excess Return ~ Market Excess Return)
for i in 1:n
    # Create a temporary DataFrame for the current asset's excess returns and the market excess return
    df = DataFrame(y = exret[:, i], x = exmkt)
    # Fit the linear regression model using GLM: y = alpha + beta * x + epsilon
    model = lm(@formula(y ~ x), df)
    # Store the coefficients (alpha and beta)
    beta[:, i] = coef(model)
end

# Calculate the mean excess return of assets
# Calculate mean along dimension 1 (columns, i.e., for each asset) and transpose to get a column vector
muExret = convert(Matrix{Float64}, mean(exret, dims=1)')  # Shape: (n, 1)

# Calculate the covariance matrix of the excess returns
Sigma = cov(exret; dims=1)  # Covariance matrix, dims=1 computes covariance between columns

# Calculate the implied market equilibrium excess returns (Pi) using CAPM: Pi = beta * Market_Excess_Return
# Use the beta coefficient (second row of the beta matrix, index 2 in Julia)
# Multiply beta by the mean of the market excess return (muExret's first element is the mean market excess return if muExret was calculated differently, but here it's mean asset excess returns).
# A more standard calculation of Pi is beta * mean(exmkt)
# Assuming the original code's intent was beta * mean(exret) element-wise, which is unusual for CAPM equilibrium.
# Let's calculate Pi as beta * mean(exmkt) as is standard in Black-Litterman.
mean_exmkt = mean(exmkt)
Pi_standard = beta[2, :] .* mean_exmkt # Element-wise multiplication of betas by mean market excess return
Pi = reshape(Pi_standard, :, 1) # Reshape to a column vector for calculations

# Uncertainty factor (Litterman and He, 1999) - often set to a small value or 1/T
tau = 1 / T

# --- Define the Investor Views (Q and P) ---
# Q is the vector of expected returns on the views, kx1 (k is the number of views)
# Example views:
Q = [0.04/12,    # View 1: Expected return for a combination of assets
     0.02/12,    # View 2: Expected return for another combination
     0.10/12]    # View 3: Expected return for a single asset
Q = reshape(Q, :, 1) # Reshape to a column vector

# P is the pick matrix, kxn (k views, n assets)
# Defines which assets are involved in each view and their weights
# Example P matrix corresponding to the example Q views:
P = [
    1 0 -1 0 0;  # View 1: Asset 1 - Asset 3
    0 -1 0 1 0;  # View 2: Asset 4 - Asset 2
    0 0 0 0 1    # View 3: Asset 5
]

# Omega: Covariance matrix of the views, kxk
# Represents the uncertainty in the investor's views. It's typically diagonal
# with variances of the view errors. A common approach is P * (tau * Sigma) * P'
Omega = P * (tau * Sigma) * P'

# --- Blend the Equilibrium returns (Pi) with the Views (Q) ---
# Calculate the Black-Litterman covariance matrix (SigmaBL)
# Formula: [(tau*Sigma)^-1 + P'*Omega^-1*P]^-1
SigmaBL = inv(inv(tau * Sigma) + P' * inv(Omega) * P)

# Calculate the Black-Litterman expected returns (muBL)
# Formula: SigmaBL * [(tau*Sigma)^-1*Pi + P'*inv(Omega)*Q]
muBL = SigmaBL * (inv(tau * Sigma) * Pi + P' * inv(Omega) * Q)

# --- Display and Save Results (Expected Returns Comparison) ---
# Create a DataFrame to compare mean excess returns, implied equilibrium returns (Pi), and Black-Litterman returns (muBL)
muTable = DataFrame(
    Asset = asset_columns, # Use the original asset names as index
    exret = round.(vec(muExret); digits=4), # Mean historical excess returns
    Pi = round.(vec(Pi); digits=4),         # Implied equilibrium returns (CAPM)
    muBL = round.(vec(muBL); digits=4)      # Black-Litterman blended returns
)
println("Comparison of Expected Returns:")
println(muTable)
# Write the table to a CSV file
XLSX.writetable("muTable.xlsx", muTable; overwrite=true) # Writing to XLSX as per other Julia scripts

# --- Implement the two efficient frontiers (Mean-Variance and Black-Litterman) ---
# Define a range of target expected returns for plotting the efficient frontiers
muR = convert(Vector{Float64}, 0:0.00001:0.015)

# Calculate the Mean-Variance efficient frontier using the mean excess returns (muExret) and covariance (Sigma)
# getuncef is a helper function assumed to be available via include("getuncef.jl")
# The getuncef function is expected to return (OptSigma, w) based on its usage here.
OptSigma, w = getuncef(muExret, Sigma, muR)

# Calculate the Black-Litterman efficient frontier using the blended returns (muBL) and adjusted covariance (Sigma + SigmaBL)
OptSigmaBL, wBL = getuncef(muBL, Sigma + SigmaBL, muR)

# --- Create Efficient Frontier Table ---
# Create a DataFrame to store the efficient frontier data (returns and standard deviations)
ef_table = DataFrame(
    r = muR,                 # Target returns
    sigma = sqrt.(OptSigma), # Standard deviation for M-V frontier
    sigmaBL = sqrt.(OptSigmaBL) # Standard deviation for B-L frontier
)

# --- Plot the Efficient Frontiers ---
# Create a plot object
p = plot(
    ef_table.sigma, ef_table.r, # Plot M-V frontier (Risk vs. Return)
    label="M-V", linewidth=1.5, color=:blue, # Added color for clarity
    xlabel="Portfolio Risk", ylabel="Portfolio Expected Return", # Axis labels
    title="Efficient Frontier", size=(800, 600), # Title and figure size
    legend=:topleft # Legend position
)
# Add the Black-Litterman efficient frontier to the plot
plot!(p, # Plot onto the existing plot p
    ef_table.sigmaBL, ef_table.r, # Plot B-L frontier (Risk vs. Return)
    label="B-L", linewidth=1.5, color=:red # Added color for clarity
)
# Set x-axis limits
xlims!(0, 0.1)

# Save the plot as an EPS file
savefig(p, "BL.png") # Specify the plot object to save

# --- Optimal Weights Calculation ---
# Define a specific target return for which to calculate optimal weights
target_return = 0.0060
# Find the index corresponding to the target return in the muR vector
# Use findfirst for the first occurrence of the target return
idx = findfirst(==(target_return), muR)

# Check if the target return was found
if idx === nothing
    println("Target return $target_return not found in muR range.")
else
    # Extract the optimal weights for the target return from the calculated weight matrices
    # w and wBL are expected to have shape (n, length(muR))
    Weights_mv = w[:, idx]' # Transpose to get a row vector (1, n)
    Weights_bl = wBL[:, idx]' # Transpose to get a row vector (1, n)

    # Vertically concatenate the weights for the two methods
    Weights = vcat(Weights_mv, Weights_bl)

    # Create a DataFrame to store the optimal weights
    Weights_df = DataFrame(Weights, :auto) # Create DataFrame from the combined weights matrix
    rename!(Weights_df, Symbol.(asset_columns)) # Rename columns to asset names
    Weights_df.Method = ["M-V", "B-L"] # Add a column to identify the method

    println("\nOptimal Weights for Target Return $(round(target_return*100, digits=2))%:")
    println(Weights_df)

    # Write the weights table to an Excel file
    XLSX.writetable("Weights_BL.xlsx", Weights_df; overwrite=true) # Specify overwrite=true
end


# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
