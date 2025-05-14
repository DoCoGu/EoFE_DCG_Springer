# MV.jl
# This script implements the Mean-Variance portfolio optimization for a set
# of risky assets and plots the efficient frontier, also calculating
# optimal portfolio weights in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, Statistics, LinearAlgebra, Plots # Added necessary packages

# Set working directory to current file location
cd(@__DIR__)

# --- Data Loading and Preparation ---

# Load asset returns from an Excel file
# Read the "Returns" sheet from the Excel file into a DataFrame
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

# Extract asset return data (assuming it starts from the 2nd column)
# Convert to a Matrix
R = Matrix(Ret[:, 2:end])
N = size(R, 2) # Number of assets

# --- Calculate Statistical Measures ---
# Calculate mean returns (expected returns) for risky assets
z = mean(R, dims=1)'  # Calculate mean along dimension 1 (columns) and transpose to get a column vector (Nx1)

# Calculate standard deviations (risk) for risky assets
sig = std(R, dims=1) # Calculate standard deviation along dimension 1 (columns)

# Calculate the covariance matrix of returns for risky assets
V = cov(R)
# Calculate the inverse of the covariance matrix
V1 = inv(V)

# --- Calculate Parameters for the Efficient Frontier ---
# Calculate the parameters A, B, C, and D from the portfolio theory
# These parameters define the shape and location of the efficient frontier
# A = z' * inv(V) * z
A = z' * V1 * z
# B = 1' * inv(V) * z
B = ones(1, N) * V1 * z # ones(1, N) creates a row vector of ones
# C = 1' * inv(V) * 1
C = ones(1, N) * V1 * ones(N) # ones(N) creates a column vector of ones
# D = A*C - B^2
D = A * C - B.^2 # Use broadcasting (.^) for element-wise squaring

# Define a range of target portfolio expected returns for plotting the efficient frontier
# Create a range from 0 to 0.015 with a step of 0.0001
mu_p = 0:0.0001:0.015

# Calculate the variance of the efficient frontier for the range of expected returns
# Formula: (1/D) * (C*mu_p^2 - 2*B*mu_p + A)
sig2_p = (1/D[]) .* (C[] .* mu_p.^2 .- 2 .* B[] .* mu_p .+ A[]) # Use broadcasting and [] to extract scalar values
# Calculate the standard deviation (risk) of the efficient frontier
sig_p = sqrt.(sig2_p) # Use broadcasting (sqrt.) for element-wise square root

# --- Plotting the Efficient Frontier ---
# Create the initial plot for the Efficient Frontier
p = plot(sig_p, mu_p, label="Efficient frontier", lw=1.5, color=:black) # Plot Risk vs. Expected Return

# Add a horizontal dashed line at the expected return of the Global Minimum Variance Portfolio (GMVP)
# GMVP Expected Return = B/C
plot!(p, [0, 1/sqrt(C[])], [B[]/C[], B[]/C[]], linestyle=:dash, label="B/C", lw=1.5, color=:black) # Plot onto plot p

# Add a vertical dashed line at the standard deviation of the Global Minimum Variance Portfolio (GMVP)
# GMVP Standard Deviation = 1/sqrt(C)
plot!(p, [1/sqrt(C[]), 1/sqrt(C[])], [0, B[]/C[]], linestyle=:dash, label="1/sqrt(C)", lw=1.5, color=:black) # Plot onto plot p

# Add scatter points for individual stocks (Risk vs. Return)
scatter!(p, sig[:], z[:], label="Stocks") # Plot onto plot p, use [:] to flatten to vectors

# Add a scatter point for the Global Minimum Variance Portfolio (GMVP)
# GMVP Standard Deviation = 1/sqrt(C), GMVP Expected Return = B/C
scatter!(p, [1/sqrt(C[])], [B[]/C[]], label="MVP", marker=:circle, ms=5) # Plot onto plot p, use [] for scalars

# Configure plot appearance
title!("Efficient Frontier") # Plot title
xlabel!("Portfolio Risk") # X-axis label
ylabel!("Portfolio Expected Return") # Y-axis label
xlims!(0, 0.1) # Set x-axis limits
plot!(p, legend=:topleft, legendfontsize=10) # Legend position and font size

# Save the plot as a PNG file
savefig(p, "MV.png") # Specify plot object and filename

# --- Optimal weights Calculation ---
# Calculate parameters for portfolio weights formula on the efficient frontier
# g and h are vectors that define the linear relationship between weights and expected return
# Formula: w_p = g + h * mu_p
g = 1/D[] * (A[] * V1 * ones(N) .- B[] * V1 * z) # Use broadcasting (.-) and [] for scalars
h = 1/D[] * (C[] * V1 * z .- B[] * V1 * ones(N)) # Use broadcasting (.-) and [] for scalars

# Define a range of target expected returns for which to calculate optimal weights
# Create a range from 0.001 to 0.05 with 10 points
mup_target = range(0.001, stop=0.05, length=10)

# Calculate optimal weights for the target expected returns
# Formula: w_p = g + h * mu_p
# This calculates weights for each target return in mup_target
wp = hcat([g .+ h .* μ for μ in mup_target]...) # Use broadcasting (.+) and (.*) and hcat to combine column vectors into a matrix

# --- Global Minimum Variance Portfolio weights ---
# Calculate the weights for the Global Minimum Variance Portfolio (GMVP)
# The GMVP has an expected return of B/C
w_mvp = g + h * (B[]/C[]) # Use [] for scalars

# --- Prepare Weights Table ---
# Define column names for the weights table
# First column is "GMVP", followed by target return percentages
col_names = ["GMVP"; string.(round.(mup_target .* 100, digits=2)) .* "%"] # Use broadcasting (.*) and string.()

# Define row names for the weights table (asset names)
row_names = names(Ret)[2:end] # Asset names from original Returns DataFrame

# Combine the GMVP weights and the optimal weights for target returns into a single matrix
all_weights = round.([w_mvp wp]; digits=3) # Vertically concatenate GMVP weights and wp matrix

# Create a DataFrame from the combined weights matrix
Weights = DataFrame(all_weights, Symbol.(col_names)) # Use Symbol.() to convert string column names to Symbols

# Add a column for the row names (Stock/Asset names)
Weights.Stock = row_names
# Reorder columns to have 'Stock' first, then the weight columns
select!(Weights, :Stock, Not(:Stock)) # Select 'Stock' first, then all other columns except the original 'Stock'

# --- Save to Excel ---
# Write the weights table to an Excel file
XLSX.writetable("Weights_MV.xlsx", Weights; overwrite=true) # Specify overwrite=true to replace if file exists


# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
