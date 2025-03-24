import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Setting
cutoff = 150
America = 1

# Import Data
Ret = pd.read_excel("../Data/IndustryPortfolios.xlsx")
Rf = pd.read_excel("../Data/FF_FactorsDaily.xlsx")
Rf = Rf['RF']
DatesReturn = pd.to_datetime(Ret.iloc[:, 0], format='%Y%m%d')

Events = pd.read_excel('../Data/Events.xlsx')

if America == 0:
    Events = Events[Events.iloc[:, 12] > cutoff].iloc[:, 3]
else:
    Events = Events[(Events.iloc[:, 12] > cutoff) & (Events["Zone"] == 1)].iloc[:, 3]

# Ret = Ret.iloc[1:, 1].values / Ret.iloc[:-1, 1].values - 1
Ret = Ret.iloc[:, 13].values

# Convert dates in serial date
Events = pd.to_datetime(Events.values.flatten(), format='%d-%m-%Y')

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

# Investment Strategy

# Buy and Hold strategy
BuyandHold = np.cumsum(Ret - Rf)

# Aviation disasters
AviationStrategy = Ret - Rf
AviationStrategy[EDummy[:, 0] == 1] = Rf[EDummy[:, 0] == 1] - Ret[EDummy[:, 0] == 1]
AviationStrategy = np.cumsum(AviationStrategy)

# Plotting
Inv_plot = plt.figure(1, figsize=(15, 10))
plt.plot(DatesReturn, BuyandHold, '-k', linewidth=1.5, label="Buy and Hold")
plt.plot(DatesReturn, AviationStrategy, '--b', linewidth=1.5, label="Aviation disasters Strategy")
plt.legend(fontsize=16, loc="best")
plt.ylabel('Cumulative return', fontsize=18)
plt.xlabel('Date', fontsize=18)
plt.tight_layout()
plt.savefig('ES_InvStrategy.eps', format='eps')
plt.show()
