using XLSX, DataFrames, Statistics, LinearAlgebra, GLM, CSV, CovarianceMatrices


function newey_west_se(X::Matrix, y::Vector, bandwidth::Int)
    n, k = size(X)
    b = X \ y
    u = y .- X * b

    S = zeros(k, k)

    for l in 0:bandwidth
        w_l = 1.0
        if l > 0
            w_l -= l / (bandwidth .+ 1)
        end

        Γ_l = zeros(k, k)
        for t in (l .+ 1):n
            xt = X[t, :]
            xt_l = X[t - l, :]
            Γ_l = Γ_l .+ xt' * xt_l * u[t] * u[t - l]
        end

        S = S .+ w_l * (Γ_l + Γ_l')
        if l == 0
            S = S .- Γ_l  # remove duplicate at lag 0
        end
    end

    XtX_inv = inv(X' * X)
    V = XtX_inv * S * XtX_inv
    return b, sqrt.(diag(V))
end


# Data Import
Ret = DataFrame(XLSX.readtable("../Data/Portfolios.xlsx", "Portfolios"))
Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

Uncertainty = DataFrame(XLSX.readtable("../Data/Uncertainty.xlsx", "Uncertainty"))
Uncertainty = convert(Vector{Float64}, Uncertainty[:, 3])
Consumption = DataFrame(XLSX.readtable("../Data/Consumption.xlsx", "Consumption"))
Consumption = convert(Vector{Float64}, Consumption[:,2])
# Compute EPU and Consumption Growth
EPU = Uncertainty[2:end,1] ./ Uncertainty[1:end-1] .- 1
C = Consumption[2:end, 1] ./ Consumption[1:end-1] .- 1

# Market returns and risk-free rate
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Excess returns
ExRet = (Matrix(Ret[2:end, 2:end])/100 .- Rf[2:end])

# Factors
Factors = hcat(Mkt, C, EPU)

# First Stage Regression
n1, n2 = size(ExRet)
nF = size(Factors, 2)

CoefAll = fill(NaN, nF, n2)
Res = fill(NaN, n1, n2)

for i in 1:n2
    df = DataFrame(Mkt=Factors[:, 1], ExRet=ExRet[:, i], C=Factors[:, 2], EPU=Factors[:, 3])
    model = lm(@formula(ExRet ~ Mkt + C + EPU), df)

    Coef = coef(model)
    CoefAll[:, i] = Coef[2:end]
    Res[:, i] = residuals(model)
end

VarCovErr = cov(Res)

# Second Stage Regression
MeanRet = convert(Matrix{Float64}, mean(ExRet, dims=1)')
Betas = convert(Matrix{Float64}, CoefAll')

df = DataFrame(hcat(Betas, MeanRet), :auto)
rename!(df, ["Mkt", "C", "EPU", "MeanRet"])

model2 = lm(@formula(MeanRet ~ 0 + Mkt + C + EPU), df)

Lambda = coef(model2)
SE = stderror(model2)
Tstat = Lambda ./ SE

# Shanken Correction
Sigma_f = cov(Factors)
BtB_inv = inv(Betas' * Betas)
VarLam = (BtB_inv * Betas' * VarCovErr * Betas * BtB_inv * (1 + Lambda' * inv(Sigma_f) * Lambda) + Sigma_f) / n1
SE_Shanken = sqrt.(diag(VarLam))
Tstat_Shanken = Lambda ./ SE_Shanken

# Time-series of cross-section regressions
LambdaFull = fill(NaN, n1, nF)
for j in 1:n1
    df = DataFrame(hcat(ExRet[j,:], Betas), :auto)
    rename!(df, ["ExRet", "Mkt", "C", "EPU"])
    model_j = lm(@formula(ExRet ~ 0 + Mkt + C + EPU), df)
    LambdaFull[j, :] = coef(model_j)'
end 

LambdaMean = mean(LambdaFull, dims=1)[:]

nF = size(LambdaFull, 2)
SE_NW = zeros(nF)

for k in 1:nF
    y = LambdaFull[:, k]
    X = ones(n1, 1)
    # Compute HAC robust variance (Newey-West with lag 2)
    b, SE = newey_west_se(X, y, 1)

    # Extract robust SE
    SE_NW[k] = SE[1]
end

Tstat_NW = LambdaMean ./ SE_NW

# Results Table
NamePort = names(Ret)[2:end]
FirstStageReg = DataFrame(Betas, :auto)
rename!(FirstStageReg, [:Mkt, :dC, :EPU])
FirstStageReg.rowid = NamePort

SecondStage = DataFrame(
    Lambda = Lambda,
    tstat = Tstat,
    t_stat_HAC = Tstat_NW,
    t_stat_Shanken = Tstat_Shanken,
    rowindex = ["Mkt", "dC", "EPU"]
)
SecondStage = stack(SecondStage, Not(:rowindex))
SecondStage = unstack(SecondStage, :variable, :value)

# Write to CSV
CSV.write("FirstStage_EPU.csv", FirstStageReg)
CSV.write("SecondStage_EPU.csv", SecondStage)
