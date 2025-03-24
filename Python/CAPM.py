import pandas as pd
import numpy as np
import statsmodels.api as sm


# Data
Ret = pd.read_excel('../Data/Returns.xlsx')
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100

ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# Testing the CAPM (One-Factor Model)

# Stocks
nS, mS = ExRet.shape
alpha = np.full((2, mS), np.nan)
beta = np.full((2, mS), np.nan)
R2 = np.full((1, mS), np.nan)
for i in range(mS):
    model = sm.OLS(ExRet[:, i], sm.add_constant(Mkt))
    results = model.fit()
    alpha[0, i] = round(results.params[0], 5)
    alpha[1, i] = round(results.tvalues[0], 5)
    beta[0, i] = round(results.params[1], 5)
    beta[1, i] = round(results.tvalues[1], 5)
    R2[0, i] = round(results.rsquared_adj, 5)

# Display Results
Stocks = pd.DataFrame(
    data= np.concatenate((alpha, beta, R2)),
    columns=Ret.iloc[:, 1:].columns,
    index=['alpha','(alpha t-stat)', 'beta',' (beta t-stat)','Adj. R2']
)

Stocks.to_csv('CAPM_Stock.csv', index=True)
