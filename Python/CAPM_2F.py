import pandas as pd
import numpy as np
import statsmodels.api as sm

# Data
# Import Stock Returns
Ret = pd.read_excel('../Data/Returns.xlsx')
Factors = pd.read_excel('../Data/FF_Factors.xlsx')
Uncertainty = pd.read_excel('../Data/Uncertainty.xlsx')


Mkt = np.array(Factors.iloc[1:, 1]) / 100
Rf = np.array(Factors.iloc[1:, 4]) / 100
ExRet = np.array(Ret.iloc[:, 1:]) - Rf.reshape(-1, 1)


VIX = np.array(Uncertainty.iloc[1:, 1].values / Uncertainty.iloc[:-1, 1].values - 1)
EPU = np.array(Uncertainty.iloc[1:, 2].values / Uncertainty.iloc[:-1, 2].values - 1)



# Testing the CAPM (One-Factor Model)
# VIX
nS, mS = ExRet.shape
alpha_VIX = np.empty((2, mS))
beta1_VIX = np.empty((2, mS))
beta2_VIX = np.empty((2, mS))
R2_VIX = np.empty((1,mS))


for i in range(mS):
    model = sm.OLS(ExRet[:, i],  sm.add_constant(np.array([Mkt,VIX]).transpose()))
    results = model.fit()
    alpha_VIX[0, i] = round(results.params[0], 5)
    alpha_VIX[1, i] = round(results.tvalues[0], 5)
    beta1_VIX[0, i] = round(results.params[1], 5)
    beta1_VIX[1, i] = round(results.tvalues[1], 5)
    beta2_VIX[0, i] = round(results.params[2], 5)
    beta2_VIX[1, i] = round(results.tvalues[2], 5)
    R2_VIX[0, i] = round(results.rsquared_adj, 5)

# Testing the CAPM (One-Factor Model)
# EPU
nS, mS = ExRet.shape
alpha_EPU = np.empty((2, mS))
beta1_EPU = np.empty((2, mS))
beta2_EPU = np.empty((2, mS))
R2_EPU = np.empty((1,mS))


for i in range(mS):
    model = sm.OLS(ExRet[:, i], sm.add_constant(np.array([Mkt,EPU]).transpose()))
    results = model.fit()
    alpha_EPU[0, i] = round(results.params[0], 5)
    alpha_EPU[1, i] = round(results.tvalues[0], 5)
    beta1_EPU[0, i] = round(results.params[1], 5)
    beta1_EPU[1, i] = round(results.tvalues[1], 5)
    beta2_EPU[0, i] = round(results.params[2], 5)
    beta2_EPU[1, i] = round(results.tvalues[2], 5)
    R2_EPU[0, i] = round(results.rsquared_adj, 5)


# Display Results
Stocks_VIX = pd.DataFrame(
    data= np.concatenate((alpha_VIX, beta1_VIX, beta2_VIX, R2_VIX)),
    index=['alpha', '(alpha t-stat)', 'beta_mkt', 'beta_mkt (t-stat)',
           'beta_vix', '(beta_vix t-stat)', 'Adj. R2'],
    columns=Ret.iloc[:, 1:].columns
)


Stocks_VIX.to_csv('CAPM_3F.csv')

Stocks_EPU = pd.DataFrame(
    data= np.concatenate((alpha_EPU, beta1_EPU, beta2_EPU, R2_EPU)),
    index=['alpha', '(alpha t-stat)', 'beta_mkt', 'beta_mkt (t-stat)',
           'beta_epu', '(beta_epu t-stat)', 'Adj. R2'],
    columns=Ret.iloc[:, 1:].columns
)


Stocks_EPU.to_csv('CAPM_3F.csv')

