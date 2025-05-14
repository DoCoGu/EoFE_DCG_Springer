# MV_rf.jl
# This script implements Mean-Variance portfolio optimization with a
# risk-free asset, deriving and plotting the Capital Market Line (CML),
# and calculates optimal portfolio weights in Julia.
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

# Load Fama-French Factors from an Excel file
# Read the "Factors" sheet from the Excel file into a DataFrame
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Extract asset return data (assuming it starts from the 2nd column)
# Convert to a Matrix
R = Matrix(Ret[:, 2:end])
N = size(R, 2) # Number of assets

# Extract risk-free rate (assuming it's in the 5th column, starting from the 2nd row)
# Convert to Float64 Vector and divide by 1200 to get monthly decimal rate
Rf = convert(Vector{Float64}, Factors[2:end, 5]) ./ 1200

# --- Calculate Statistics ---
# Calculate mean returns (expected returns) for risky assets
z = mean(R, dims=1)'  # Calculate mean along dimension 1 (columns) and transpose to get a column vector (Nx1)

# Calculate standard deviations (risk) for risky assets
sig = std(R, dims=1)

# Calculate the covariance matrix of returns for risky assets
V = cov(R)
# Calculate the inverse of the covariance matrix
V1 = inv(V)

# Calculate the mean of the risk-free rate
mean_rf = mean(Rf)

# Calculate excess expected returns for risky assets (Expected Return - Risk-Free Rate)
excess_z = z .- mean_rf

# Calculate H, which is related to the squared Sharpe Ratio of the tangency portfolio
# H = (z - rf*1)' * inv(S) * (z - rf*1)
H = (excess_z)' * V1 * excess_z
sqrtH = sqrt(H[]) # Extract the scalar value from the 1x1 matrix H and take the square root

# --- Efficient Frontier with Risk-Free Asset (Capital Market Line) ---
# Define a range of target portfolio expected returns for plotting the CML
# Start from the mean risk-free rate up to 0.009, with 100 points
mu_p = range(mean_rf, stop=0.009, length=100)

# Calculate the variance of portfolios on the Capital Market Line
# Formula: (1/H) * (mu_p - rf)^2
sig2_p = (1/H[]) .* (mu_p .- mean_rf).^2 # Use broadcasting (.*) for element-wise operations
# Calculate the standard deviation (risk) of portfolios on the CML
sig_p = sqrt.(sig2_p) # Use broadcasting (sqrt.) for element-wise square root

# --- Calculate Tangency Point ---
# Calculate the expected return of the Tangency Portfolio
# mu_t = rf + H / (1' * inv(S) * (z - rf*1))
mu_t = mean_rf + H[] / (ones(1, N) * V1 * excess_z)[] # Use [] to extract scalar from 1x1 matrices

# Calculate the standard deviation of the Tangency Portfolio
# sig_t = sqrt(H) / (1' * inv(S) * (z - rf*1))
sig_t = sqrtH / (ones(1, N) * V1 * excess_z)[] # Use [] to extract scalar from 1x1 matrices

# --- Efficient Frontier without Risk-Free Asset (for comparison) ---
# Recalculate parameters A, B, C, D for the risky asset efficient frontier
A = z' * V1 * z
B = z' * V1 * ones(N) # ones(N) creates a column vector of ones
C = ones(1, N) * V1 * ones(N) # ones(1, N) creates a row vector of ones
D = A*C - B.^2 # Use broadcasting (.^) for element-wise squaring

# Define a range of target expected returns for the risky asset efficient frontier
mu_p2 = range(0.001, stop=0.008, length=100)

# Calculate the variance of the risky asset efficient frontier
# Formula: (1/D) * (C*mu_p2^2 - 2*B*mu_p2 + A)
sig2_pp = (1/D[]) .* (C[] .* mu_p2.^2 .- 2 .* B[] .* mu_p2 .+ A[]) # Use broadcasting and [] for scalars
# Calculate the standard deviation (risk) of the risky asset efficient frontier
sig_pp = sqrt.(sig2_pp) # Use broadcasting (sqrt.)

# Set standard deviations below a certain threshold to NaN for cleaner plotting
# This is often done to show only the upper (efficient) part of the frontier
sig_pp[mu_p2 .< 0.005] .= NaN # Use broadcasting (.=) to assign NaN to elements meeting the condition

# --- Plotting ---
# Create the initial plot for the Capital Market Line (CML)
p = plot(sig_p, mean_rf .+ sqrtH .* sig_p, lw=1.5, label="Efficient Frontier (CML)", color=:blue) # Use broadcasting (.+) and (.*)

# Add the risky asset efficient frontier to the plot
plot!(p, sig_pp, mu_p2, lw=1.5, label="Risky Asset Frontier", color=:black) # Plot onto plot p, added label

# Add a scatter point for the Tangency Portfolio
scatter!(p, [sig_t], [mu_t], label="Tangency Portfolio", marker=:circle, ms=5, color=:red) # Plot onto plot p, added label

# Configure plot appearance
plot!(p, xlabel="Portfolio Risk", ylabel="Portfolio Expected Return", title="Capital Market Line", # Axis labels and title
      legend=:topleft, legendfontsize=12, xlims=(0, 0.05)) # Legend position, font size, and x-axis limits

# Add a horizontal line at the expected return of the Tangency Portfolio
plot!(p, [0, sig_t], [mu_t, mu_t], linestyle=:dash, color=:black, lw=1.5, label="") # Plot onto plot p, added label="" to exclude from legend

# Calculate the lower part of the CML (borrowing/lending at risk-free rate)
Ep_ = mean_rf .- sqrtH .* sig_p # Use broadcasting (.-) and (.*)
# Set values below a certain threshold to NaN for cleaner plotting (to show only the line below rf)
Ep_[Ep_ .< mean_rf - 0.002] .= NaN # Adjusted threshold slightly for clarity

# Add the lower part of the CML to the plot
plot!(p, sig_p, Ep_, lw=1.5, color=:blue, label="") # Plot onto plot p, added label=""

# Add a vertical line from the minimum of the lower CML to the Tangency Portfolio's expected return
# Need to find the minimum of the non-NaN values in Ep_
min_Ep_val = minimum(filter(!isnan, Ep_))
plot!(p, [sig_t, sig_t], [min_Ep_val, mu_t], linestyle=:dash, color=:black, lw=1.5, label="") # Plot onto plot p, added label=""


# Save the plot as a PNG file
savefig(p, "MV_rf.png") # Specify plot object and filename

# --- Optimal weights Calculation ---
# Define a range of target expected returns for which to calculate optimal weights on the CML
mup_target = range(0.001, stop=0.05, length=10)

# Calculate optimal weights for portfolios on the CML (combination of risk-free and tangency)
# Formula: w_p = inv(S) * (z - rf*1) * (mu_p - rf) / H
# This calculates weights for each target return in mup_target
wp = hcat([(V1 * excess_z) * (μ - mean_rf) / H[] for μ in mup_target]...) # Use hcat to combine column vectors into a matrix

# Calculate the weights for the Tangency Portfolio (risky assets only)
# Formula: w_T = inv(S) * (z - rf*1) / (1' * inv(S) * (z - rf*1))
wT = (V1 * excess_z) / (ones(1, N) * V1 * excess_z)[] # Use [] to extract scalar

# Calculate the weight allocated to the risk-free asset (1 - sum of risky asset weights)
# sum(hcat(wT, wp), dims=1) calculates the sum across rows (dimension 1) for each column (portfolio)
rf_weights = 1 .- sum(hcat(wT, wp), dims=1) # Use broadcasting (.-)

# --- Prepare Weights Table ---
# Define column names for the weights table
col_names = ["wT"; string.(round.(mup_target .* 100, digits=2)) .* "%"] # Column names: Tangency weights, then target returns in percentage

# Define row names for the weights table (asset names + Rf)
row_names = [names(Ret)[2:end]; "Rf"] # Asset names from original Returns DataFrame, plus "Rf"

# Combine all calculated weights into a single matrix
# wT and wp are weights for risky assets. rf_weights are weights for the risk-free asset.
all_weights = round.([wT wp; rf_weights]; digits=3) # Vertically concatenate risky asset weights and RF weights

# Create a DataFrame from the combined weights matrix
Weights = DataFrame(all_weights, Symbol.(col_names)) # Use Symbol.() to convert string column names to Symbols

# Add a column for the row names (Stock/Asset names and Rf)
Weights.Stock = row_names
# Reorder columns to have 'Stock' first, then the weight columns
select!(Weights, :Stock, Not(:Stock)) # Select 'Stock' first, then all other columns except the original 'Stock'

# --- Save to Excel ---
# Write the weights table to an Excel file
XLSX.writetable("Weights_MVrf.xlsx", Weights; overwrite=true) # Specify overwrite=true to replace if file exists


# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
