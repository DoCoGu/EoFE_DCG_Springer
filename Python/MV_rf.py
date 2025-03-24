import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import data
Ret = pd.read_excel('Returns.xlsx')
Factors = pd.read_excel('FF_Factors.xlsx')

Rf = np.array(Factors.iloc[1:, 4]) / 1200

R = np.array(Ret.iloc[:, 1:])
N = R.shape[1]

z = np.mean(R, axis=0)
sig = np.std(R, axis=0)
V = np.cov(R, rowvar=False)
V1 = np.linalg.inv(V)

H = np.dot((z - np.mean(Rf)), np.dot(V1, (z - np.mean(Rf))))
sqrtH = np.sqrt(H)

mu_p = np.linspace(np.mean(Rf), 0.009, 100)

A = np.dot(z.T, np.dot(V1, z))
B = np.dot(z.T, np.dot(V1, np.ones((N, 1))))
C = np.dot(np.ones((1, N)), np.dot(V1, np.ones((N, 1))))
D = A * C - B ** 2

# Variance and standard deviation
sig2_p = 1 / H * (mu_p - np.mean(Rf)) ** 2
sig_p = np.sqrt(sig2_p)

mu_t = np.mean(Rf) + H / (np.ones((1, N)) @ V1 @ (z - np.ones((1, N)) * np.mean(Rf)).T)
sig_t = sqrtH / (np.ones((1, N)) @ V1 @ (z - np.ones((1, N)) * np.mean(Rf)).T)

A = np.dot(z.T, np.dot(V1, z))
B = np.dot(z.T, np.dot(V1, np.ones((N, 1))))
C = np.dot(np.ones((1, N)), np.dot(V1, np.ones((N, 1))))
D = A * C - B ** 2
mu_p = np.linspace(0.001, 0.008, 100)

# Variance and standard deviation
sig2_p = 1 / D * (C * mu_p ** 2 - 2 * B * mu_p + A)
sig_pp = np.sqrt(sig2_p)
sig_pp[0,mu_p < 0.005] = np.nan
x = np.vstack((sig_pp, mu_p))



# Plot
fig = plt.figure(2)
plt.plot(sig_p, np.mean(Rf) + sqrtH * sig_p, '-b', linewidth=1.5)
plt.plot(sig_pp.transpose(), mu_p, '-k', linewidth=1.5)
plt.scatter(sig_t, mu_t, marker='o', color='k')
plt.title('Capital market line', fontsize=16)
plt.xlabel('Portfolio Risk', fontsize=16)
plt.ylabel('Portfolio Expected Return', fontsize=16)
plt.legend(["Efficient frontier", "Portfolio", "Tangency port.", "Stocks", "Autoupdate", "off"], fontsize=16, loc="upper left")
plt.xlim([0, 0.05])
plt.ylim([0.002, 0.009])
plt.axhline(mu_t, linestyle='--', color='k', linewidth=1.5)
plt.axvline(sig_t, linestyle='--', color='k', linewidth=1.5)
plt.savefig("MV_rf.eps", format="eps")
plt.show()



mup = np.linspace(0.001, 0.05, 10)
wp = V1 @ ((z - np.ones((1, N))).T * np.mean(Rf)) * (mup - np.mean(Rf)) / H

# Tangency portfolio weights
wT = (V1 @ ((z - np.ones((1, N))).T * np.mean(Rf))) / (np.ones((N, 1)).T @ V1 @ ((z - np.ones((1, N))).T * np.mean(Rf)))

# Weight on the risk-free asset
rf = 1 - np.sum(np.concatenate((wT, wp),axis=1), axis=0)

weight_data = np.vstack([np.concatenate((wT, wp),axis=1), rf])
weight_data = np.round(weight_data, 3)
columns = ["wT"] + [f"{round(mup_val * 100, 2)}%" for mup_val in mup]
index = list(Ret.columns[1:]) + ["Rf"]
Weights = pd.DataFrame(weight_data, columns=columns, index=index)

Weights.to_excel("Weights_MVrf.xlsx", index=True, header=True)
