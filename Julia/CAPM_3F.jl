using XLSX, DataFrames, GLM, Statistics

# === Load Data ===

Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))


# === Data Preparation ===

# Extract Mkt and Rf
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Excess returns
ExRet = Matrix(Ret[:, 2:end]) .- Rf[2:end]

Factors = Matrix(Factors[2:end, 2:4]) ./ 100  # size: T x 3

# === Run Regressions ===

nS, mS = size(ExRet)

# Initialize result matrices

nS, mS = size(ExRet)
alpha = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_1 = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_2 = Array{Union{Missing, Float64}}(undef, 2, mS)
beta_3 = Array{Union{Missing, Float64}}(undef, 2, mS)
R2 = Array{Union{Missing, Float64}}(undef, mS)

for i in 1:mS
    df = DataFrame(Mkt=Factors[:, 1], ExRet=ExRet[:, i], SMB=Factors[:, 2], HML=Factors[:, 3])
    model = lm(@formula(ExRet ~ Mkt + SMB + HML), df)

    coefs = coef(model)
    ctable = coeftable(model)
    tstat = ctable.cols[3]

    alpha[1, i] = coefs[1]
    alpha[2, i] = tstat[1]
    beta_1[1, i] = coefs[2]
    beta_1[2, i] = tstat[2]
    beta_2[1, i] = coefs[3]
    beta_2[2, i] = tstat[3]
    beta_3[1, i] = coefs[4]
    beta_3[2, i] = tstat[4]

    R2[i] = adjr2(model)
end


result_matrix = round.([alpha; beta_1; beta_2; beta_3; R2'], digits=5)

# === Format Results for Excel ===

col_names = names(Ret)[2:end]
row_names = [
    "alpha", "(alpha t-stat)",
    "beta_mkt", "beta_mkt (t-stat)",
    "beta_smb", "(beta_smb t-stat)",
    "beta_hml", "(beta_hml t-stat)",
    "Adj. R2"
]


Results = DataFrame(result_matrix', :auto)
rename!(Results, Symbol.(row_names))
Results = permutedims(Results)  # Transpose for correct orientation
rename!(Results, col_names)

rename!(Results, Symbol.(col_names))
Results = hcat(DataFrame(Statistic = row_names), Results)

XLSX.writetable("CAPM_3F.xlsx", Results; overwrite=true)