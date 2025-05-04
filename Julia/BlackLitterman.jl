using XLSX
using DataFrames
using LinearAlgebra
using Statistics
using GLM
using Plots

include("getuncef.jl")  

# The Black-Litterman model - Application

# Import and prepare data
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

Ret[!, "Mkt-RF"] = Factors[2:end, "Mkt-RF"] ./ 100
Ret[!, "rf"] = Factors[2:end, "RF"] ./ 100

# Calculate excess returns
asset_columns = names(Ret)[2:end-2]
exret = Matrix(Ret[:, asset_columns]) .- Ret.rf
exmkt = Ret[!, "Mkt-RF"]

# Compute the equilibrium returns
T, n = size(exret)
beta = zeros(2, n)

for i in 1:n
    model = lm(@formula(y ~ x), DataFrame(y = exret[:, i], x = exmkt))
    beta[:, i] = coef(model)
end

muExret = convert(Matrix{Float64},mean(exret, dims=1)')  # 5×1 vector
Sigma = cov(exret; dims=1)      # Covariance matrix
Pi = beta[2, :] .* muExret
Pi = reshape(Pi, :, 1)          # Reshape for calculations
tau = 1 / T

# Define the Q and P
Q = [0.04/12, 0.02/12, 0.10/12]
Q = reshape(Q, :, 1)

P = [
    1 0 -1 0 0;
    0 -1 0 1 0;
    0 0 0 0 1
]

# Omega: Covariance matrix of the views
Omega = P * (tau * Sigma) * P'

# Blend the equilibrium returns Pi with the views
SigmaBL = inv(inv(tau * Sigma) + P' * inv(Omega) * P)
muBL = SigmaBL * (inv(tau * Sigma) * Pi + P' * inv(Omega) * Q)

# Table of muExret, Pi, muBL
muTable = DataFrame(
    Asset = asset_columns,
    exret = round.(vec(muExret); digits=4),
    Pi = round.(vec(Pi); digits=4),
    muBL = round.(vec(muBL); digits=4)
)
XLSX.writetable("muTable.csv", muTable; overwrite=true)


# Implement the two efficient frontiers
muR =  convert(Vector{Float64}, 0:0.00001:0.015)

OptSigma, w = getuncef(muExret, Sigma, muR)
OptSigmaBL, wBL = getuncef(muBL, Sigma + SigmaBL, muR)

# Create efficient frontier table
ef_table = DataFrame(
    r = muR,
    sigma = sqrt.(OptSigma),
    sigmaBL = sqrt.(OptSigmaBL)
)

# Plot the efficient frontiers
plot(
    ef_table.sigma, ef_table.r,
    label="M-V", linewidth=1.5,
    xlabel="Portfolio Risk", ylabel="Portfolio Expected Return",
    title="Efficient Frontier", size=(800, 600)
)
plot!(
    ef_table.sigmaBL, ef_table.r,
    label="B-L", linewidth=1.5
)
xlims!(0, 0.1)
savefig("BL.eps")

# Optimal weights
idx = findfirst(==(0.0060), muR)
Weights = vcat(
    w[:, idx]',
    wBL[:, idx]'
)
Weights_df = DataFrame(Weights, :auto)
rename!(Weights_df, Symbol.(asset_columns))
Weights_df.Method = ["M-V", "B-L"]

XLSX.writetable("Weights_BL.xlsx", Weights_df)