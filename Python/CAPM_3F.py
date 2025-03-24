import pandas as pd
import numpy as np
import statsmodels.api as sm


# Data
Ret = pd.read_excel('../Data/Returns.xlsx')
Factors = pd.read_excel('../Data/FF_Factors.xlsx')

FF = np.array(Factors.iloc[1:, 1:4]) / 100

Rf = np.array(Factors.iloc[1:, 4]) / 100

ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)

# Testing the CAPM (3-Factor Model)

# Stocks
nS, mS = ExRet.shape
alpha = np.empty((2, mS))
beta1 = np.empty((2, mS))
beta2 = np.empty((2, mS))
beta3 = np.empty((2, mS))
R2 = np.empty((1,mS))


for i in range(mS):
    model = sm.OLS(ExRet[:, i], sm.add_constant(FF))
    results = model.fit()
    alpha[0, i] = round(results.params[0], 5)
    alpha[1, i] = round(results.tvalues[0], 5)
    beta1[0, i] = round(results.params[1], 5)
    beta1[1, i] = round(results.tvalues[1], 5)
    beta2[0, i] = round(results.params[2], 5)
    beta2[1, i] = round(results.tvalues[2], 5)
    beta3[0, i] = round(results.params[3], 5)
    beta3[1, i] = round(results.tvalues[3], 5)
    R2[0, i] = round(results.rsquared_adj, 5)




# Display Results

Stocks = pd.DataFrame(
    data= np.concatenate((alpha, beta1, beta2, beta3, R2)),
    index=['alpha', '(alpha t-stat)', 'beta_mkt', 'beta_mkt (t-stat)',
           'beta_smb', '(beta_smb t-stat)', 'beta_hml', '(beta_hml t-stat)', 'Adj. R2'],
    columns=Ret.iloc[:, 1:].columns
)


Stocks.to_csv('CAPM_3F.csv')
