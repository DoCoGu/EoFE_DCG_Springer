using XLSX, DataFrames, Statistics, LinearAlgebra, Plots

# Load returns data
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))
R = Matrix(Ret[:, 2:end])
N = size(R, 2)

# Statistical measures
z = mean(R, dims=1)'  # expected returns
sig = std(R, dims=1)  # std deviation
V = cov(R)
V1 = inv(V)

A = z' * V1 * z
B = z' * V1 * ones(N)
C = ones(1, N) * V1 * ones(N)
D = A*C - B.^2

mu_p = 0:0.0001:0.015

# Portfolio variance and standard deviation
sig2_p = (1/D) .* (C .* mu_p.^2 .- 2 .* B .* mu_p .+ A)
sig_p = sqrt.(sig2_p)

# Plot Efficient Frontier
plot(sig_p, mu_p, label="Efficient frontier", lw=1.5, color=:black)
plot!([0, 1/sqrt(C[])], [B[]/C[], B[]/C[]], linestyle=:dash, label="B/C", lw=1.5, color=:black)
plot!([1/sqrt(C[]), 1/sqrt(C[])], [0, B[]/C[]], linestyle=:dash, label="1/sqrt(C)", lw=1.5, color=:black)
scatter!(sig[:], z[:], label="Stocks")
scatter!([1/sqrt(C[])], [B/C], label="MVP", marker=:circle, ms=5)
title!("Efficient Frontier")
xlabel!("Portfolio Risk")
ylabel!("Portfolio Expected Return")
xlims!(0, 0.1)
plot!(legend=:topleft, fontsize=10)
savefig("MV.png")

# Optimal weights
g = 1/D[] * (A[] * V1 * ones(N) .- B[] * V1 * z)
h = 1/D[] * (C[] * V1 * z .- B[] * V1 * ones(N))

mup = range(0.001, stop=0.05, length=10)
wp = hcat([g .+ h .* μ for μ in mup]...)

# Global Minimum Variance Portfolio weights
w_mvp = g + h * (B/C)

# Prepare weights table
col_names = ["GMVP"; string.(round.(mup .* 100, digits=2)) .* "%"]
row_names = names(Ret)[2:end]


all_weights = round.([w_mvp wp]; digits=3)

Weights = DataFrame(all_weights, Symbol.(col_names))

Weights.rowval = row_names
select!(rename!(Weights, :rowval => :Stock), :Stock, Not(:Stock))

# Write to Excel
XLSX.writetable("Weights_MV.xlsx", Weights; overwrite=true)

