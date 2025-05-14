# getuncef.jl
# This function calculates the unconstrained efficient frontier for a given
# set of asset expected returns and covariance matrix in Julia.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

using XLSX, DataFrames, Statistics, LinearAlgebra, Plots # Added necessary packages (though Plots is not used in this function itself)

"""
    getuncef(mu::Matrix, S::Matrix, muR)

Calculates the unconstrained efficient frontier for a given set of asset
expected returns and covariance matrix.

Parameters:
    mu : Matrix
        Expected returns (column vector, shape (k, 1)).
    S : Matrix
        Covariance matrix (shape (k, k)).
    muR : Float64 or Vector
        Target returns (can be a single value or array of values).

Returns:
    OptSigma : Float64 or Vector
        Optimal portfolio variance corresponding to target returns.
    w : Matrix
        Weights for the efficient frontier portfolios (shape (k, length(muR)) if muR is array).
    # Note: Original MATLAB function also returned gmv and tradeoff, but this Julia
    # version's return signature only includes OptSigma and w based on the provided code.
"""
function getuncef(mu::Matrix, S::Matrix, muR)

    k = size(S, 1) # Get the number of assets
    inv_S = inv(S) # Calculate the inverse of the covariance matrix

    # Calculate the parameters A, B, C, and D from the portfolio theory
    # These parameters define the shape and location of the efficient frontier
    # A = mu' * inv(S) * mu
    A = (mu' * inv_S * mu)[1] # Extract the scalar value from the 1x1 matrix
    # B = 1' * inv(S) * mu
    B = (ones(k)' * inv_S * mu)[1] # ones(k)' creates a row vector of ones, extract scalar
    # C = 1' * inv(S) * 1
    C = (ones(k)' * inv_S * ones(k))[1] # ones(k) creates a column vector of ones, extract scalar
    # D = A*C - B^2
    D = A * C - B^2

    # Calculate the variance of the efficient frontier for given expected returns (muR)
    # Formula: (1/D) * (C*muR^2 - 2*B*muR + A)
    OptSigma = (C .* muR.^2 .- 2 .* B .* muR .+ A) ./ D # Use broadcasting for element-wise operations

    # Calculate parameters related to the tradeoff between risk and return
    # These are components used in calculating efficient frontier weights
    # d = A * inv(S) * 1 - B * inv(S) * mu
    d = A * (inv_S * ones(k)) - B * (inv_S * mu)
    # E = C * inv(S) * mu - B * inv(S) * 1
    E = C * (inv_S * mu) - B * (inv_S * ones(k))

    # Calculate the weights of the efficient frontier portfolios for given expected returns (muR)
    # Formula: w = d/D + E/D * muR
    # Use broadcasting for element-wise operations with muR (which can be a scalar or vector)
    w = d ./ D .+ (E ./ D) .* muR' # Use muR' to handle both scalar and vector muR correctly

    # Return the calculated values
    # Note: The return signature matches the provided function code (OptSigma, w)
    return OptSigma, w
end

# This Julia script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
