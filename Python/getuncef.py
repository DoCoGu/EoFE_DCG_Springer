# getuncef.py
# This function calculates the unconstrained efficient frontier for a given
# set of asset expected returns and covariance matrix in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

import numpy as np

def getuncef(mu, S, muR):
    """
    Python Code
    Function: getuncef - get unconstrained efficient frontier
    
    Parameters:
        mu : ndarray
            Expected returns (column vector).
        S : ndarray
            Covariance matrix.
        muR : float or ndarray
            Target returns (can be a single value or array of values).
            
    Returns:
        muR : float or ndarray
            Target returns (same as input).
        OptSigma : float or ndarray
            Optimal portfolio variance corresponding to target returns.
        gmv : ndarray
            Global minimum variance portfolio [expected return, variance].
        tradeoff : ndarray
            Tradeoff metrics [risk-return ratio, risk per unit of variance].
        w : ndarray
            Weights for the efficient frontier portfolios.
    """
    k = S.shape[1]
    
    inv_S = np.linalg.inv(S)
    
    A = mu.T @ inv_S @ mu
    B = np.ones((k, 1)).T @ inv_S @ mu
    C = np.ones((k, 1)).T @ inv_S @ np.ones((k, 1))
    D = A * C - B**2
    
    OptSigma = (C * muR**2 - 2 * B * muR + A) / D
    
    
    d = A * inv_S @ np.ones((k, 1)) - B * inv_S @ mu
    E = C * inv_S @ mu - B * inv_S @ np.ones((k, 1))
    w = d / D + (E / D) * muR
    
    return OptSigma, w

# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
