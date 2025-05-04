using XLSX, DataFrames, Statistics, GLM, Plots, Dates

# -- Data --

# Import Stock Returns
Ret = DataFrame(XLSX.readtable("../Data/Returns.xlsx", "Returns"))

Factors = DataFrame(XLSX.readtable("../Data/FF_Factors.xlsx", "Factors"))


# Extract Mkt and Rf
Mkt = convert(Vector{Float64}, Factors[2:end, 2]) ./ 100
Rf = convert(Vector{Float64}, Factors[:, 5]) ./ 100

# Excess returns
ExRet = Matrix(Ret[:, 2:end]) .- Rf[2:end]

# -- Regress --

# Rolling Window
rw = 40
nS, mS = size(ExRet)
alpha = fill(NaN, nS - rw + 1, mS)
beta = fill(NaN, nS - rw + 1, mS)

for t in 1:(nS - rw + 1)
    for i in 1:mS
        df = DataFrame(ExRet=ExRet[t:t+rw-1, i], Mkt=Mkt[t:t+rw-1])
        model = lm(@formula(ExRet ~ Mkt), df)
        alpha[t, i] = coef(model)[1]
        beta[t, i] = coef(model)[2]
    end
end

# -- Plots --

# Generate plot dates
dates = Ret[rw:end, 1]  # Assuming dates are in the first column
dates = Date.(dates)
day_range = LinRange(0, Dates.value(dates[end] - dates[1]), 20)
Dates_for_plot = [start_date + Day(round(Int, d)) for d in day_range]

date_labels = Dates.format.(DateTime.(Dates_for_plot), dateformat"u-yyyy")

# Alpha plot
plot(1:nS - rw + 1, alpha[:, 1], label=names(Ret)[2], lw=1.5, color=:red)
for i in 2:mS
    plot!(1:nS - rw + 1, alpha[:, i], label=names(Ret)[i+1], lw=1.5)
end
plot!(1:nS - rw + 1, zeros(nS - rw + 1), ls=:dash, color=:black, label="")
xticks!(range(1, stop=nS - rw + 1, length=20), date_labels)
plot!(xtickfont=font(8))
plot!(xrotation=45)
title!("Alpha")
savefig("alpha_S")

# Beta plot
plot(1:nS - rw + 1, beta[:, 1], label=names(Ret)[2], lw=1.5, color=:red)
for i in 2:mS
    plot!(1:nS - rw + 1, beta[:, i], label=names(Ret)[i+1], lw=1.5)
end
plot!(1:nS - rw + 1, ones(nS - rw + 1), ls=:dash, color=:black, label="")
xticks!(range(1, stop=nS - rw + 1, length=20), date_labels)
plot!(xtickfont=font(8))
plot!(xrotation=45)
title!("Beta")
savefig("Beta_S")
