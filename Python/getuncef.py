import numpy as np

def getuncef(mu, S, muR):
    """
    ESGSF2021 - michele.costola@unive.it
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