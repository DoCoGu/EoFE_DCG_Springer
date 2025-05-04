using XLSX, DataFrames, Statistics, LinearAlgebra, Plots

# Load data from Excel
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

R = Matrix(Ret[:, 2:end])
N = size(R, 2)

Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 1200

# --- Stats ---
z = mean(R, dims=1)'        # Expected returns (Nx1)
sig = std(R, dims=1)
V = cov(R)
V1 = inv(V)

mean_rf = mean(Rf)
excess_z = z .- mean_rf
H = (excess_z)' * V1 * excess_z
sqrtH = sqrt(H[])

# Efficient frontier with risk-free asset (Capital Market Line)
mu_p = range(mean_rf, stop=0.009, length=100)
sig2_p = (1/H[]) .* (mu_p .- mean_rf).^2
sig_p = sqrt.(sig2_p)

# Tangency point
mu_t = mean_rf + H[] / (ones(1, N) * V1 * excess_z)[]
sig_t = sqrtH / (ones(1, N) * V1 * excess_z)[]

# Efficient frontier without risk-free asset
A = z' * V1 * z
B = z' * V1 * ones(N)
C = ones(1, N) * V1 * ones(N)
D = A*C - B.^2

mu_p2 = range(0.001, stop=0.008, length=100)
sig2_pp = (1/D[]) .* (C[] .* mu_p2.^2 .- 2 .* B[] .* mu_p2 .+ A[])
sig_pp = sqrt.(sig2_pp)
sig_pp[mu_p2 .< 0.005] .= NaN

# --- Plotting ---
plot(sig_p, mean_rf .+ sqrtH .* sig_p, lw=1.5, label="Efficient Frontier (CML)", color=:blue)
plot!(sig_pp, mu_p2, lw=1.5, label="Portfolio", color=:black)
scatter!([sig_t], [mu_t], label="Tangency Portfolio", marker=:circle, ms=5, color=:red)
plot!(xlabel="Portfolio Risk", ylabel="Portfolio Expected Return", title="Capital Market Line",
      legend=:topleft, legendfontsize=12, xlims=(0, 0.05))

# Horizontal line at mu_t
plot!([0, sig_t], [mu_t, mu_t], linestyle=:dash, color=:black, lw=1.5, label="")
# Vertical line at sig_t
plot!([sig_t, sig_t], [minimum(filter(!isnan, Ep_)), mu_t], linestyle=:dash, color=:black, lw=1.5, label="")

# Plot expected excess returns line
Ep_ = mean_rf .- sqrtH .* sig_p
Ep_[Ep_ .< -0.002] .= NaN
plot!(sig_p, Ep_, lw=1.5, color=:blue, label="")

# Vertical line from Ep_ min to mu_t

savefig("MV_rf.png")

# --- Optimal weights ---
mup = range(0.001, stop=0.05, length=10)
wp = hcat([(V1 * excess_z) * (μ - mean_rf) / H[] for μ in mup]...)

# Tangency portfolio weights
wT = (V1 * excess_z) / (ones(N)' * V1 * excess_z)[]

# Weight on risk-free asset
rf_weights = 1 .- sum(hcat(wT, wp), dims=1)


# --- Prepare Weights Table ---
col_names = ["wT"; string.(round.(mup .* 100, digits=2)) .* "%"]
row_names = [names(Ret)[2:end]; "Rf"]
all_weights = round.([wT wp; rf_weights]; digits=3)

Weights = DataFrame(all_weights, Symbol.(col_names))
Weights.rowval = row_names
select!(rename!(Weights, :rowval => :Stock), :Stock, Not(:Stock))


# --- Save to Excel ---
XLSX.writetable("Weights_MVrf.xlsx", Weights; overwrite=true)

