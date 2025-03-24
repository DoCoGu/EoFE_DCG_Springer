import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from getuncef import getuncef

# The Black-Litterman model - Application

# Import and prepare data
ret = pd.read_excel('../Data/Returns.xlsx')  # 5 assets
ff = pd.read_excel('../Data/FF_Factors.xlsx')  # FamaFrench factors

ret['Mkt-RF'] = ff['Mkt-RF'][1:].values / 100
ret['rf'] = ff['RF'][1:].values / 100

# Calculate excess returns
exret = ret.iloc[:, 1:-2].values - ret['rf'].values[:,None]  # Excess returns
exmkt = ret['Mkt-RF'].values  # Excess market returns

# Compute the equilibrium returns
T, n = exret.shape
beta = np.empty((2, n))

for i in range(n):
    mdl = sm.OLS(exret[:, i], sm.add_constant(exmkt)).fit()
    beta[:, i] = mdl.params

muExret = exret.mean(axis=0)
Sigma = np.cov(exret, rowvar=False)  # Covariance matrix of the excess returns
Pi = beta[1, :] * muExret  # Market equilibrium returns
Pi = Pi[:, None]  # Reshape for calculations
tau = 1 / T  # Uncertainty factor

# Define the Q and P
Q = np.array([0.04 / 12, 0.02 / 12, 0.10 / 12])[:, None]  # Expected returns on views (monthly)
P = np.array([
    [1, 0, -1, 0, 0],  # Relative view: Asset 1 vs Asset 3
    [0, -1, 0, 1, 0],  # Relative view: Asset 2 vs Asset 4
    [0, 0, 0, 0, 1]    # Absolute view: Asset 5
])

# Omega: Covariance matrix of the views
Omega = P @ (tau * Sigma) @ P.T

# Blend the equilibrium returns Pi with the views
SigmaBL = np.linalg.inv(np.linalg.inv(tau * Sigma) + P.T @ np.linalg.inv(Omega) @ P)
muBL = SigmaBL @ (np.linalg.inv(tau * Sigma) @ Pi + P.T @ np.linalg.inv(Omega) @ Q)

# Table of muExret, Pi, muBL
muTable = pd.DataFrame(
    np.round(np.hstack([muExret[:, None], Pi, muBL]), 4),
    columns=['exret', 'Pi', 'muBL'],
    index=ret.columns[1:-2]
)
muTable.to_csv('muTable.csv', index_label='Asset')

# Implement the two efficient frontiers
muR = np.arange(0, 0.015, 0.00001)



OptSigma, w = getuncef(muExret[:,None], Sigma, muR)
OptSigmaBL, wBL = getuncef(muBL, Sigma+SigmaBL, muR)

# Create efficient frontier table
ef_table = pd.DataFrame(
    {
        'r': muR.flatten(),
        'sigma': np.sqrt(OptSigma.flatten()),
        'sigmaBL': np.sqrt(OptSigmaBL.flatten())
    }
)

# Plot the efficient frontiers
plt.figure(figsize=(10, 6))
plt.plot(ef_table['sigma'], ef_table['r'], '-b', label='M-V', linewidth=1.5)
plt.plot(ef_table['sigmaBL'], ef_table['r'], '-r', label='B-L', linewidth=1.5)
plt.title('Efficient Frontier', fontsize=16)
plt.xlim([0, 0.1])
plt.xlabel('Portfolio Risk', fontsize=16)
plt.ylabel('Portfolio Expected Return', fontsize=16)
plt.legend(loc='upper left', fontsize=14)
plt.savefig('BL.eps', format='eps')
plt.show()

# Optimal weights
idx = np.where(muR == 0.0060)[0][0]
Weights = pd.DataFrame(
    np.vstack([w[:, idx], wBL[:, idx]]),
    columns=ret.columns[:5],
    index=['M-V', 'B-L']
)
Weights.to_excel('Weights_BL.xlsx', index_label='Method')
