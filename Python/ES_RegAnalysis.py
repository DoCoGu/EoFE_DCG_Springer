import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import statsmodels.api as sm

# Setting
maxlag = 3
cutoff = 150
fD = 3
cons = 1
America = 1

# Import Data
Ret = pd.read_excel("../Data/IndustryPortfolios.xlsx")
DatesReturn = pd.to_datetime(Ret.iloc[:, 0], format='%Y%m%d')

Events = pd.read_excel('../Data/Events.xlsx')

if America == 0:
    Events = Events[Events.iloc[:, 12] > cutoff].iloc[:, 3]
else:
    Events = Events[(Events.iloc[:, 12] > cutoff) & (Events["Zone"] == 1)].iloc[:, 3]

Ret = Ret.iloc[:, 13].values

# Convert dates in serial date
Events = pd.to_datetime(Events.values.flatten(), format='%d-%m-%Y')

# Preallocate vectors of weekday dummies
dow = np.full((len(DatesReturn), 4), np.nan)

# Create dummy for weekdays
for d in range(4):
    for i in range(len(DatesReturn)):
        dow[i, d] = 1 if DatesReturn[i].weekday() == d else 0

# Create events dummy
Match = np.isin(DatesReturn, Events)

EDummy = np.zeros((len(DatesReturn)+3, 3))

for i in range(len(DatesReturn)):
    if Match[i]:
        if i + 1 < len(DatesReturn):
            EDummy[i + 1, 0] = 1
        if i + 2 < len(DatesReturn):
            EDummy[i + 2, 1] = 1
        if i + 3 < len(DatesReturn):
            EDummy[i + 3, 2] = 1

EDummy = EDummy[3:, :]

T = np.zeros(len(Ret))
dday = DatesReturn.dt.day
dmonth = DatesReturn.dt.month

T[(dmonth == 1) & (dday >= 1) & (dday <= 5)] = 1

# Prepare regression data
RetLag = np.column_stack([np.roll(Ret, lag) for lag in range(1, maxlag + 1)])
RetLag[:maxlag, :] = np.nan

X = np.column_stack([RetLag, dow, T, EDummy[:, :fD]])

# Adding a constant for the intercept in statsmodels
X = sm.add_constant(X[maxlag:])  # Add intercept
Ret_trimmed = Ret[maxlag:]

# Fit the linear regression model
mdl = sm.OLS(Ret_trimmed, X).fit()

# Extract coefficients and standard errors
coeff = np.round(mdl.params, 5)
se = mdl.bse  # Standard errors

Sig = coeff / se
res = np.vstack([coeff, Sig])

varnames = ["Mon", "Tue", "Wed", "Thu", "TaxDays"]
rString = [f"R_t-{lag}" for lag in range(1, maxlag + 1)]
eString = [f"Event +{i}" for i in range(1, fD + 1)]
var = ["Constant"] + rString + varnames + eString

# Create results table
tab = pd.DataFrame(res, columns=var, index=["Coeff", "Sig"])

if America == 0:
    tab.to_excel('RegressionAnalysis_All.xlsx')
else:
    tab.to_excel('RegressionAnalysis_America.xlsx')