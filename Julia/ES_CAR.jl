using XLSX, CSV, DataFrames, Dates, Statistics, LinearAlgebra, GLM, Plots

# Settings
America = 0      # Use only aviation disasters in America
cutoff = 150     # Cutoff for number of casualties
start_CAR = -1
end_CAR = 5


start_CAPM = -250
end_CAPM = -50

# Import Data
Portfolios = DataFrame(XLSX.readtable("../Data/IndustryPortfolios.xlsx", "Portfolios"))
Factors = DataFrame(XLSX.readtable("../Data/FF_FactorsDaily.xlsx", "Factors"))
Events = DataFrame(XLSX.readtable("../Data/Events.xlsx", "Events"))

# Filter Events
if America == 0
    Events = Events[(Events[:, 13] .> cutoff), 4]
else
    Events = Events[(Events[:, 13] .> cutoff) .& (Events[:, "Zone"] .== 1), 4]
end

# Convert returns
MktRet = Factors[:, 2] ./ 100
Rf = Factors[:,5] ./ 100
Ret = Portfolios[:, 14] ./ 100
ExRet = Ret .- Rf
ExRet = reshape(ExRet, :, 1)  # Make it 2D

# Convert dates
Events_dates = Date.(Events[:, 1], dateformat"dd-mm-yyyy")
DatesReturn = Date.(string.(Factors[:, 1]), dateformat"yyyymmdd")

n, nS = size(ExRet)

# Match events to return dates
Dates_vec = zeros(Int, length(Events_dates))
for i in 1:length(Events_dates)
    event_date = Events_dates[i]
    idx = findfirst(x -> x == event_date, DatesReturn)
    while isnothing(idx)
        event_date += Day(1)
        idx = findfirst(x -> x == event_date, DatesReturn)
    end
    Dates_vec[i] = idx
end

Dates_event = Dates_vec[Dates_vec .> abs(start_CAPM)]
Dates_event = unique(Dates_event)
nE = length(Dates_event)

# Event Study: Estimate alpha and beta
alpha = fill(NaN, nE, nS)
beta = fill(NaN, nE, nS)

for i in 1:nE
    for k in 1:nS
        X = MktRet[Dates_event[i] + start_CAPM : Dates_event[i] + end_CAPM]
        Y = ExRet[Dates_event[i] + start_CAPM : Dates_event[i] + end_CAPM, k]
        df = DataFrame(X = X, Y = Y)
        model = lm(@formula(Y ~ X), df)
        alpha[i, k] = coef(model)[1]
        beta[i, k] = coef(model)[2]
    end
end

# Predicted returns
len_CAR = abs(start_CAR) + abs(end_CAR) + 1
PredRet = fill(NaN, len_CAR, nE, nS)

for t in 1:len_CAR
    for i in 1:nE
        for k in 1:nS
            PredRet[t, i, k] = alpha[i, k] + beta[i, k] * MktRet[Dates_event[i] + start_CAR - 1 + t]
        end
    end
end

# Observed returns
ObsRet_agg = fill(NaN, len_CAR, nE, nS)
for i in 1:nE
    for t in 1:len_CAR
        for k in 1:nS
            ObsRet_agg[t, i, k] = ExRet[Dates_event[i] + start_CAR - 1 + t, k]
        end
    end
end

# Abnormal and Cumulative Returns
AbnRet = ObsRet_agg .- PredRet
CAR = mapslices(x -> cumsum(skipmissing(x)), AbnRet, dims=1)
CAAR = mean(CAR, dims=2)[:]
CAAR .*= 100

# Plot
date = collect(start_CAR:end_CAR)
zero_vec = zeros(length(date))

plot(date, CAAR, lw=2, label="CAR", color=:blue)
plot!(date, zero, lw=2, linestyle=:dash, color=:black, label="")
vline!([0], lw=2, color=:black, linestyle=:solid, label="")
xlabel!("Days Relative to Events")
ylabel!("CAR")
xlims!(start_CAR, end_CAR)
ylims!(minimum(CAAR) - 0.001, maximum(CAAR) + 0.001)
title!("Cumulative Abnormal Returns (CAR)")
plot!(legend=:topleft)
savefig(America == 0 ? "CAAR_All" : "CAAR_America")
