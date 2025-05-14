# getuncef.R
# This function calculates the unconstrained efficient frontier for a given
# set of asset expected returns and covariance matrix in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Note: This function typically requires libraries like 'stats' (for cov, mean, solve)
# and 'base' (for matrix operations, rep, t, etc.), which are usually loaded by default in R.

#' getuncef - Get Unconstrained Efficient Frontier
#'
#' Calculates the unconstrained efficient frontier for a given set of asset
#' expected returns and covariance matrix.
#'
#' @param mu A vector or matrix (k x 1) of expected returns.
#' @param S A covariance matrix (k x k).
#' @param muR A numeric vector or scalar of target returns for the efficient frontier.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{muR}: The target returns (same as input).
#'     \item \code{OptSigma}: Optimal portfolio variance corresponding to target returns.
#'     \item \code{w}: Weights for the efficient frontier portfolios (matrix, k x length(muR)).
#'   }
#'   Note: Original MATLAB function also returned gmv and tradeoff, but this R
#'   version's return signature only includes muR, OptSigma, and w based on the provided code.
#'
#' @examples
#' # Example usage (requires hypothetical data)
#' # mu_example <- c(0.01, 0.015, 0.02) # Expected returns
#' # S_example <- matrix(c(0.005, 0.001, 0.002, 0.001, 0.008, 0.003, 0.002, 0.003, 0.01), 3, 3) # Covariance matrix
#' # muR_example <- seq(0.01, 0.025, length.out = 50) # Target returns
#' # ef_results <- getuncef(mu_example, S_example, muR_example)
#' # plot(sqrt(ef_results$OptSigma), ef_results$muR, type = 'l',
#' #      xlab = 'Portfolio Standard Deviation', ylab = 'Portfolio Expected Return')
#'
getuncef <- function(mu, S, muR) {
  # Ensure mu is treated as a column vector/matrix if it's a simple vector
  if (!is.matrix(mu)) {
    mu <- as.matrix(mu)
  }

  k <- ncol(S)  # Number of assets (from the covariance matrix dimensions)

  # Precompute the inverse of the covariance matrix for efficiency
  inv_S <- as.matrix(solve(S))

  # Create a column vector of ones with the same number of rows as assets
  ones_k <- rep(1, k)

  # Calculate the parameters A, B, C, and D from the portfolio theory
  # These parameters define the shape and location of the efficient frontier
  # A = mu' * inv(S) * mu
  A <- as.vector(t(mu) %*% inv_S %*% mu) # Use as.vector to extract the scalar value
  # B = 1' * inv(S) * mu
  B <- as.vector(t(ones_k) %*% inv_S %*% mu) # Use as.vector to extract the scalar value
  # C = 1' * inv(S) * 1
  C <- as.vector(t(ones_k) %*% inv_S %*% ones_k) # Use as.vector to extract the scalar value
  # D = A*C - B^2
  D <- as.vector(A * C - B^2) # Use as.vector to ensure D is a scalar

  # Calculate the variance of the efficient frontier for given expected returns (muR)
  # Formula: (1/D) * (C*muR^2 - 2*B*muR + A)
  # Use broadcasting for element-wise operations (muR^2, 2*B*muR, etc.)
  OptSigma <- (C * muR^2 - 2 * B * muR + A) / D

  # Calculate parameters related to the tradeoff between risk and return
  # These are components used in calculating efficient frontier weights
  # d = A * inv(S) * 1 - B * inv(S) * mu
  d <- A * inv_S %*% ones_k - B * inv_S %*% mu
  # E = C * inv(S) * mu - B * inv(S) * 1
  E <- C * inv_S %*% mu - B * inv_S %*% ones_k

  # Calculate the weights of the efficient frontier portfolios for given expected returns (muR)
  # Formula: w = d/D + E/D * muR
  # d/D and E/D are vectors (k x 1). muR is a vector (1 x length(muR)) if input is vector, or scalar.
  # outer() is used to perform the element-wise multiplication of E/D (vector) by muR (vector or scalar).
  # drop() removes redundant dimensions if muR was a scalar.
  # as.vector(d/D) ensures d/D is a vector before addition.
  w <- drop(as.vector(d/D) + outer(as.vector(E) / D, muR)) # Corrected outer multiplication and use as.vector for E/D

  # Return the calculated results in a list
  # The list contains the target returns, optimal variances, and portfolio weights.
  return(list(muR = muR, OptSigma = OptSigma, w = w))
}


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
