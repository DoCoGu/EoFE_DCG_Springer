import pandas as pd
import numpy as np
import statsmodels.api as sm


# Import Stock Returns
Ret = pd.read_excel('../Data/PortfoliosLong.xlsx')
FFactors = pd.read_excel('../Data/FF_FactorsLong.xlsx')
Temperature = pd.read_excel('../Data/Temperature.xlsx')
Consumption = pd.read_excel('../Data/ConsumptionLong.xlsx')

T = np.diff(np.log(Temperature.iloc[:, 1]))
C = Consumption.iloc[1:, 1].values / Consumption.iloc[:-1, 1].values - 1

Mkt = FFactors.iloc[1:, 1].values / 100
Rf = FFactors.iloc[1:, 4].values / 100


ExRet = Ret.iloc[1:, 1:].values/100 - Rf.reshape(-1,1)

Factors = pd.DataFrame({'Mkt': Mkt, 'C': C, 'T': T})


# First-stage regression
n1, n2 = ExRet.shape
nF = Factors.shape[1]

CoefAll = np.empty((nF, n2))
Res = np.empty((n1, n2))


for i in range(n2):
    model = sm.OLS(ExRet[:, i], sm.add_constant(Factors))
    results = model.fit()
    CoefAll[:,i] = results.params[1:]
    Res[:, i] = results.resid


VarCovErr = np.cov(Res.T)

# Second-stage regression
MeanRet = np.mean(ExRet, axis=0)

# Betas
Betas = CoefAll.T

# Cross-sectional regression
model = sm.OLS(MeanRet, Betas).fit()
SE = model.bse
Lambda = model.params
Tstat = Lambda / SE

# Shanken correction
Sigma_f = np.cov(Factors, rowvar=False)

B = Betas
BtB_inv = np.linalg.inv(B.T @ B)
correction = 1 + Lambda.T @ np.linalg.inv(Sigma_f) @ Lambda
VarLam = BtB_inv @ B.T @ VarCovErr @ B @ BtB_inv * correction + Sigma_f
VarLam = VarLam / n1

SE_Shanken = np.sqrt(np.diag(VarLam))
Tstat_Shanken = Lambda / SE_Shanken

# Time-series of cross-sectional regressions
nF = Betas.shape[1]
LambdaFull = np.full((n1, nF), np.nan)

for j in range(n1):
    MeanRet_j = ExRet[j, :]
    model_j = sm.OLS(MeanRet_j, Betas).fit()
    LambdaFull[j, :] = model_j.params

LambdaMean = np.mean(LambdaFull, axis=0)

# HAC-corrected standard errors (Newey-West)
X = np.ones((n1, 1))
hac_cov = np.zeros(nF)

for k in range(nF):
    y = LambdaFull[:, k]
    model_k = sm.OLS(y, X).fit(cov_type='HAC', cov_kwds={'maxlags': 2})
    hac_cov[k] = model_k.cov_params()

SE_NW = np.sqrt(hac_cov)
Tstat_NW = LambdaMean / SE_NW


## Results
NamePort = Ret.columns[1:]
FirstStageReg = pd.DataFrame(CoefAll.T, columns=['Mkt', 'dC', 'T'], index=NamePort)
SecondStage = pd.DataFrame(np.array([Lambda, Tstat, Tstat_NW, Tstat_Shanken]).T, columns=['Lambda','t-stat', 't-stat HAC', 't-stat Shanken'], index=['Mkt', 'dC', 'EPU']).T

FirstStageReg.to_csv('FirstStage_T.csv')
SecondStage.to_csv('SecondStage_T.csv')


