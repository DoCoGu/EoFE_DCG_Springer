getuncef <- function(mu, S, muR) {
  k <- ncol(S)  # Number of assets
  
  # Precompute some values for efficiency
  inv_S <- as.matrix(solve(S))
  ones_k <- rep(1, k)
  
  A <-as.vector(t(mu) %*% inv_S %*% mu)
  B <- as.vector(t(ones_k) %*% inv_S %*% mu)
  C <- as.vector(t(ones_k) %*% inv_S %*% ones_k)
  D <- as.vector(A * C - B^2)
  
  # Efficient frontier
  OptSigma <- (C * muR^2 - 2 * B * muR + A)/D


  # Weights for each return on the efficient frontier
  d <- A * inv_S %*% ones_k - B * inv_S %*% mu
  E <- C * inv_S %*% mu - B * inv_S %*% ones_k
  
  w <- drop(as.vector(d/D) + outer((as.vector(E) / D), t(muR)))

  
  
  
  # Return all the results
  return(list(muR = muR, OptSigma = OptSigma, w = w))
}