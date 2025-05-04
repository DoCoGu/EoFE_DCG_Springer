using XLSX, DataFrames, GLM, Statistics

# === Data Import ===

# Load data from Excel
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

Uncertainty = DataFrame(XLSX.readtable("../Data/Uncertainty.xlsx", "Uncertainty"))



# === Data Preparation ===

# Convert VIX and EPU to percentage changes
VIX = Vector(Uncertainty[2:end, 2] ./ Uncertainty[1:end-1,2] .-1 )

EPU = Vector(Uncertainty[2:end, 3] ./ Uncertainty[1:end-1,3] .-1 )



# Extract Mkt and Rf
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Excess returns
ExRet = Matrix(Ret[:, 2:end]) .- Rf[2:end]

# === Model Estimation ===

nS, mS = size(ExRet)
alpha_VIX = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_1_VIX = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_2_VIX = Array{Union{Missing, Float64}}(undef, 2, mS)
R2_VIX = Array{Union{Missing, Float64}}(undef, mS)


for i in 1:mS
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i], VIX = VIX)
    model = lm(@formula(ExRet ~ Mkt + VIX), df)

    coefs = coef(model)
    ctable = coeftable(model)
    tstat = ctable.cols[3]  # 4th column = t values



    alpha_VIX[1, i] = coefs[1]
    alpha_VIX[2, i] = tstat[1]
    beta_1_VIX[1, i] = coefs[2]
    beta_1_VIX[2, i] = tstat[2]
    beta_2_VIX[1, i] = coefs[3]
    beta_2_VIX[2, i] = tstat[3]

    R2_VIX[i] = adjr2(model)
end

alpha_EPU = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_1_EPU = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_2_EPU = Array{Union{Missing, Float64}}(undef, 2, mS)
R2_EPU = Array{Union{Missing, Float64}}(undef, mS)


for i in 1:mS
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i], EPU = EPU)
    model = lm(@formula(ExRet ~ Mkt + EPU), df)

    coefs = coef(model)
    ctable = coeftable(model)
    tstat = ctable.cols[3]  # 4th column = t values



    alpha_EPU[1, i] = coefs[1]
    alpha_EPU[2, i] = tstat[1]
    beta_1_EPU[1, i] = coefs[2]
    beta_1_EPU[2, i] = tstat[2]
    beta_2_EPU[1, i] = coefs[3]
    beta_2_EPU[2, i] = tstat[3]

    R2_EPU[i] = adjr2(model)
end

result_VIX = round.([alpha_S; beta_1; beta_2; R2_S'], digits=5)

row_names = ["alpha", "(alpha t-stat)", "MKt", "(MKt t-stat)","VIX", "(VIX t-stat)", "Adj. R2"]
col_names = names(Ret)[2:end]

Stocks_VIX = DataFrame(result_VIX', :auto)
rename!(Stocks_VIX, Symbol.(row_names))
Stocks_VIX = permutedims(Stocks_VIX)  # Transpose for correct orientation
rename!(Stocks_VIX, col_names)

rename!(Stocks_VIX, Symbol.(col_names))
Stocks_VIX = hcat(DataFrame(Statistic = row_names), Stocks_VIX)

result_EPU = round.([alpha_S; beta_1; beta_2; R2_S'], digits=5)

row_names = ["alpha", "(alpha t-stat)", "MKt", "(MKt t-stat)","EPU", "(EPU t-stat)", "Adj. R2"]
col_names = names(Ret)[2:end]

Stocks_EPU = DataFrame(result_EPU', :auto)
rename!(Stocks_EPU, Symbol.(row_names))
Stocks_EPU = permutedims(Stocks_EPU)  # Transpose for correct orientation
rename!(Stocks_EPU, col_names)

rename!(Stocks_EPU, Symbol.(col_names))
Stocks_EPU = hcat(DataFrame(Statistic = row_names), Stocks_EPU)


XLSX.writetable("CAPM_VIX.xlsx", Stocks_VIX; overwrite=true)
XLSX.writetable("CAPM_EPU.xlsx", Stocks_EPU; overwrite=true)
