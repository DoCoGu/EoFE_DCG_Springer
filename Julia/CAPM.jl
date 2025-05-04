using DataFrames, CSV, GLM, XLSX, Statistics

# === Data ===

# Import Stock Returns
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))


Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))

# Extract Market Return and Risk-Free Rate, and convert from percentages
Mkt = Factors[2:end, 2] ./ 100  # Market return
Rf = Factors[2:end, 5] ./ 100  # Risk-free rate

# Excess Returns
ExRet = Matrix(Ret[:, 2:end]) .- Rf

# === Testing the CAPM (One-Factor Model) ===

nS, mS = size(ExRet)
alpha_S = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_S = Array{Union{Missing, Float64}}(undef, 2, mS)
R2_S = Array{Union{Missing, Float64}}(undef, mS)

for i in 1:mS
    df = DataFrame(Mkt=Mkt, ExRet=ExRet[:, i])
    model = lm(@formula(ExRet ~ Mkt), df)

    coefs = coef(model)
    ctable = coeftable(model)
    tstat = ctable.cols[3]  # 4th column = t values



    alpha_S[1, i] = coefs[1]
    alpha_S[2, i] = tstat[1]
    beta_S[1, i] = coefs[2]
    beta_S[2, i] = tstat[2]
    R2_S[i] = adjr2(model)
end

# === Display Results ===

# Prepare results table
result_matrix = round.([alpha_S; beta_S; R2_S'], digits=5)

row_names = ["alpha", "(alpha t-stat)", "beta", "(beta t-stat)", "Adj. R2"]
col_names = names(Ret)[2:end]

Results = DataFrame(result_matrix', :auto)
rename!(Results, Symbol.(row_names))
Results = permutedims(Results)  # Transpose for correct orientation
rename!(Results, col_names)

rename!(Results, Symbol.(col_names))
Results = hcat(DataFrame(Statistic = row_names), Results)

XLSX.writetable("CAPM_Stock.xlsx", Results; overwrite=true)