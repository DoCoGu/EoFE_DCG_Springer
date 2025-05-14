# ES_InvestmentStrategy.py
# This script analyzes different investment strategies, including buy-and-hold
# and an event-driven strategy based on aviation disasters, in Python.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025


using CSV, DataFrames, XLSX, Dates, Statistics, Plots

# Set the working directory to current file location
cd(@__DIR__)

# Setting
cutoff = 150
America = 0
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


# Buy and Hold strategy: Cumulative excess return over time
# Calculate excess returns (Asset Return - Risk-Free Rate)
BuyandHold = cumsum(Ret .- Rf)

# Aviation disasters event-driven strategy
# This strategy replaces the excess return on the day after an aviation disaster
# (where EDummy[:, 0] == 1) to be (Rf - Ret) instead of (Ret - Rf).
# This implies selling the risky asset and holding the risk-free asset on event days + 1.
# Create a copy of the excess returns to modify for the strategy
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

# This Python script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895