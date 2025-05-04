function getuncef(mu::Matrix, S::Matrix, muR)
    """
    Function: getuncef - get unconstrained efficient frontier

    Parameters:
        mu : Vector
            Expected returns (column vector).
        S : Matrix
            Covariance matrix.
        muR : Float64 or Vector
            Target returns (can be a single value or array of values).

    Returns:
        OptSigma : Float64 or Vector
            Optimal portfolio variance corresponding to target returns.
        w : Matrix
            Weights for the efficient frontier portfolios.
    """

    k = size(S, 1)
    inv_S = inv(S)

    A = (mu' * inv_S * mu)[1]
    B = (ones(k)' * inv_S * mu)[1]
    C = (ones(k)' * inv_S * ones(k))[1]
    D = A * C - B^2

    OptSigma = (C .* muR.^2 .- 2 .* B .* muR .+ A) ./ D

    d = A * (inv_S * ones(k)) - B * (inv_S * mu)
    E = C * (inv_S * mu) - B * (inv_S * ones(k))
    w = d ./ D .+ (E ./ D) .* muR'

    return OptSigma, w
end