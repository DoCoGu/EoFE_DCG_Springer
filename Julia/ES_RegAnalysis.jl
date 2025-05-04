using XLSX, DataFrames, Dates, GLM, StatsModels, ShiftedArrays

# Settings
maxlag = 3
cutoff = 150
fD = 3
America = 1

# --- Import Data ---

# Read returns and dates
Portfolios = DataFrame(XLSX.readtable("../Data/IndustryPortfolios.xlsx", "Portfolios"))
DatesReturn = Date.(string.(Portfolios[:, 1]), dateformat"yyyymmdd")

# Read events
Events = DataFrame(XLSX.readtable("../Data/Events.xlsx", "Events"))

# Filter events
if America == 0
    Events = Events[Events[:, 13] .> cutoff, 4]
else
    Events = Events[(Events[:, 13] .> cutoff) .& (Events[!, "Zone"] .== 1), 4]
end

Events = Date.(Events[:, 1], dateformat"dd-mm-yyyy")

# Get return series (Transport industry column 14)
Ret = Portfolios[:, 14]
nObs = size(Ret, 1)

# --- Weekday Dummies (Mon to Thu) ---
dow = zeros(Int, nObs, 4)
for d = 1:4
    for i = 1:nObs
        if dayofweek(DatesReturn[i]) == d  # weekday 2=Mon, ..., 5=Thu
            dow[i, d] = 1
        end
    end
end

# --- Event Dummies ---
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
# --- Tax Day Dummy (Jan 1–5) ---
T = zeros(Int, nObs)
dday = day.(DatesReturn)
dmonth = month.(DatesReturn)
T[(dmonth .== 1) .& (dday .>= 1) .& (dday .<= 5)] .= 1

# --- Lagged Returns ---
RetLag = hcat([lag(Ret[:, 1], l) for l in 1:maxlag]...)
RetLag[1:maxlag, :] .= 0  # Replace initial NaNs with zeros

# --- Design Matrix ---
X = convert(Matrix{Float64},hcat(RetLag, dow, T, EDummy[:, 1:fD]))
varnames = vcat(["Rt$i" for i in 1:maxlag], ["Mon", "Tue", "Wed", "Thu", "TaxDays"], ["Event$i" for i in 1:fD])
df = DataFrame(X, Symbol.(varnames))
df.Returns = convert(Vector{Float64},Ret)
# --- Regression ---
nvars = size(X, 2)


# 2. Create DataFrame with named predictors
formula = Term(:Returns) ~ sum(Term.(Symbol.(varnames)))
model = lm(formula, df[4:end,:])

# --- Results ---
coeff = round.(coef(model), digits=5)
se = stderror(model)
Sig = coeff ./ se
res = vcat(coeff', Sig')

# --- Variable Names ---
rStrings = ["R_{t-$(i)}" for i in 1:maxlag]
varnames = ["Mon", "Tue", "Wed", "Thu", "TaxDays"]
eStrings = ["Event +$(i)" for i in 1:fD]
all_vars = ["Constant"; rStrings; varnames; eStrings]

# --- Output Table ---
results_df = DataFrame(reshape(res, 2, :), :auto)
rename!(results_df, Symbol.(all_vars))
results_df[!, :Row] = ["Coefficient", "t-Stat"]
results_df = select(results_df, :Row, Not(:Row))

# --- Save Output ---
output_file = America == 0 ? "RegressionAnalysis_All.xlsx" : "RegressionAnalysis_America.xlsx"
XLSX.writetable(output_file, "Results" => results_df)
