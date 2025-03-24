import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt


# Data
Ret = pd.read_excel('../Data/Returns.xlsx')
Factors = pd.read_excel('../Data/FF_Factors.xlsx')
Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100

ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# Testing the CAPM (One-Factor Model)
# Rolling-Window Size
rw = 40

nS, mS = ExRet.shape
alpha = np.empty((nS - rw + 1, mS))
beta = np.empty((nS - rw + 1, mS))

for t in range (nS - rw + 1):
    for i in range(mS):
        model = sm.OLS(ExRet[t:t+rw-1, i], sm.add_constant(Mkt[t:t+rw-1]))
        results = model.fit()
        alpha[t, i] = round(results.params[0], 5)
        beta[t, i] = round(results.params[1], 5)

## Plots
Dates_for_plot = pd.date_range(start=Ret.iloc[rw-1, 0], end=Ret.iloc[-1, 0], periods=20)
Dates_for_plot = [d.strftime('%b-%Y') for d in Dates_for_plot]

# alpha
alpha_plot = plt.figure(1)
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 0], 'r', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 1], 'k', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 2], 'b', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 3], 'g', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), alpha[:, 4], 'm', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), np.zeros(nS - rw + 1), '--k')
plt.xticks(np.linspace(1, nS - rw + 1, 20), Dates_for_plot, rotation=45, fontsize=14)
plt.xlim([1, nS - rw  + 1])
plt.legend(Ret.columns[1:], loc='best')
plt.title('Alpha')
plt.savefig('alpha_S.png', dpi=300)

# beta
beta_plot = plt.figure(2)
plt.plot(np.arange(1, nS - rw + 2), beta[:, 0], 'r', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), beta[:, 1], 'k', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), beta[:, 2], 'b', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), beta[:, 3], 'g', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), beta[:, 4], 'm', linewidth=1.5)
plt.plot(np.arange(1, nS - rw + 2), np.ones(nS - rw + 1), '--k')
plt.xticks(np.linspace(1, nS - rw + 1, 20), Dates_for_plot, rotation=45, fontsize=14)
plt.xlim([1, nS - rw])
plt.legend(Ret.columns[1:], loc='best')
plt.title('Beta')
plt.savefig('beta_S.png', dpi=300)
