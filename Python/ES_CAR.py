import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Settings
America = 1  # Use only aviation disasters in America
cutoff = 150  # Cutoff for number of casualties
start_CAR = -1
end_CAR = 5

start_CAPM = -250
end_CAPM = -50

# Import Data
Portfolios = pd.read_excel("../Data/IndustryPortfolios.xlsx")
Factors = pd.read_excel("../Data/FF_FactorsDaily.xlsx")
Events = pd.read_excel("../Data/Events.xlsx")

if America == 0:
    Events = Events[Events.iloc[:, 12] > cutoff].iloc[:, 3]
else:
    Events = Events[(Events.iloc[:, 12] > cutoff) & (Events.iloc[:, -1] == 1)].iloc[:, 3]



# Convert Factors to arrays
MktRet = np.array(Factors.iloc[:, 1]) / 100
Rf = np.array(Factors.iloc[:, 4]) / 100

# Convert Portfolios to arrays
Ret = np.array(Portfolios.iloc[:, 13]) / 100  # Transport industry
ExRet = Ret - Rf

# Convert dates to datetime objects
Events = [datetime.strptime(str(date), '%d-%m-%Y') for date in Events]
DatesReturn = [datetime.strptime(str(date), '%Y%m%d') for date in Factors.iloc[:, 0]]

n = ExRet.shape

Dates = []
for event_date in Events:
    while event_date not in DatesReturn:
        event_date += timedelta(days=1)
    Dates.append(DatesReturn.index(event_date))

Dates = np.array(Dates)
Dates = Dates[Dates > abs(start_CAPM)]
Dates = np.unique(Dates)
nE = len(Dates)


# Preallocate arrays for CAPM coefficients
alpha = np.empty((nE, 1))
alpha[:] = np.nan
beta = np.empty((nE, 1))
beta[:] = np.nan

# Compute the CAPM model for each stock and around each event
for i in range(nE):
    start_index = Dates[i] + start_CAPM
    end_index = Dates[i] + end_CAPM

    # Fit the linear regression model
    X = MktRet[start_index:end_index+1]
    y = ExRet[start_index:end_index+1]
    X = sm.add_constant(X)  # Add a constant term to the independent variable
    model = sm.OLS(y, X)
    results = model.fit()

    # Store the estimated coefficients
    alpha[i, 0] = results.params[0]
    beta[i, 0] = results.params[1]
    
        
# Preallocate array for CAPM predicted returns
PredRet = np.empty((abs(start_CAR) + abs(end_CAR) + 1, nE))
PredRet[:] = np.nan
# Get predicted returns for each stock and around each event, for the CAR period
for t in range(abs(start_CAR) + abs(end_CAR)+1):
    for i in range(nE):
        PredRet[t, i] = alpha[i, 0] + beta[i, 0] * MktRet[Dates[i] + start_CAR  + t]

# Preallocate array for observed returns
ObsRet_agg = np.empty((abs(start_CAR) + abs(end_CAR) + 1, nE))
ObsRet_agg[:] = np.nan

# Get observed returns for each stock and around each event, for the CAR period
for i in range(nE):
    for t in range(abs(start_CAR) + abs(end_CAR) + 1):
        ObsRet_agg[t, i] = ExRet[Dates[i] + start_CAR + t]

# Get abnormal returns (Observed ret. - Predicted ret.)
AbnRet = ObsRet_agg - PredRet

# Get cumulative abnormal returns
CAR = np.nancumsum(AbnRet, axis=0)

# Calculate CAAR (average of CAR)
CAAR = np.nanmean(CAR, axis=1) * 100
CAAR = np.squeeze(CAAR)


# Plot

# Get dates vector for plot
date = np.arange(start_CAR, end_CAR + 1).reshape(-1, 1)

zero = np.zeros((abs(start_CAR) + abs(end_CAR) + 1, 1))

CAAR_Plot = plt.figure(1)
plt.plot(date, CAAR, '-b', date, zero, '--k', linewidth=2)
plt.axvline(x=0, color='black', linestyle='-', linewidth=2)
plt.ylabel('CAR', fontsize=18)
plt.xlabel('Days Relative to Events', fontsize=18)
plt.xlim(start_CAR, end_CAR)
plt.ylim(min(CAAR) - 0.001, max(CAAR) + 0.001)
plt.legend(['CAR'])
plt.xticks(rotation=90)
plt.rcParams["figure.figsize"] = [10, 10]

if America == 0:
    plt.savefig('CAAR_All.eps', format='eps')
else:
    plt.savefig('CAAR_America.eps', format='eps')



