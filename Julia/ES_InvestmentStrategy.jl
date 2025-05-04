using CSV, DataFrames, XLSX, Dates, Statistics, Plots

# Setting
cutoff = 150
America = 1
fD = 3

# Import Data
Portfolios = DataFrame(XLSX.readtable("../Data/IndustryPortfolios.xlsx", "Portfolios"))
Factors = DataFrame(XLSX.readtable("../Data/FF_FactorsDaily.xlsx", "Factors"))
Events = DataFrame(XLSX.readtable("../Data/Events.xlsx", "Events"))
Rf = convert(Vector{Float64}, Factors[:, 5])
DatesReturn = Date.(string.(Portfolios[:, 1]), dateformat"yyyymmdd")


if America == 0
    Events = Events[Events[:, 13] .> cutoff, 4]
else
    Events = Events[(Events[:, 13] .> cutoff) .& (Events[!, "Zone"] .== 1), 4]
end

Events = Date.(Events[:, 1], dateformat"dd-mm-yyyy")


# Get return series (Transport industry column 14)
Ret = Portfolios[:, 14]
nObs = size(Ret, 1)

# --- Event Dummies ---
EDummy = zeros(Int, nObs+1, 1)
match_idx = findall(x -> x in Events, DatesReturn)

EDummy = zeros(Int, nObs+3, fD)
match_idx = findall(x -> x in Events, DatesReturn)

for i in match_idx
    for d = 1:fD
        if i + d <= nObs
            EDummy[i + d, d] = 1
        end
    end
end
EDummy = EDummy[4:end,:]


# Investment Strategy
BuyandHold = cumsum(Ret .- Rf)

AviationStrategy = copy(Ret .- Rf)
aviation_idx = EDummy[:, 1] .== 1
AviationStrategy[aviation_idx] .= Rf[aviation_idx] .- Ret[aviation_idx]

AviationStrategy = cumsum(AviationStrategy)

# Plotting
plot(DatesReturn, BuyandHold, label="Buy and Hold", linewidth=1.5, color=:black)
plot!(DatesReturn, AviationStrategy, label="Aviation disasters Strategy", linewidth=1.5, linestyle=:dash, color=:blue)
xlabel!("Date")
ylabel!("Cumulative return")
plot!(legend=:best, size=(1000, 600), legendfontsize=12, guidefontsize=14)

savefig("ES_InvStrategy.png")
