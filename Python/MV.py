import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Clear the workspace
plt.close('all')

# Read returns data from file
Ret = pd.read_excel('Returns.xlsx', engine='openpyxl', sheet_name=0)

R = Ret.iloc[:, 1:].values
N = R.shape[1]

z = np.mean(R, axis=0).reshape(-1, 1)
sig = np.std(R, axis=0)
V = np.cov(R, rowvar=False)
V1 = np.linalg.inv(V)

A = z.T @ V1 @ z
B = z.T @ V1 @ np.ones((N, 1))
C = np.ones((1, N)) @ V1 @ np.ones((N, 1))
D = A * C - B**2

mu_p = np.arange(0, 0.0151, 0.0001)

# Variance and standard deviation
sig2_p = (1 / D) * (C * mu_p**2 - 2 * B * mu_p + A)
sig_p = np.sqrt(sig2_p)


# Plot
# Create a figure
p = plt.figure(1)
p.canvas.manager.full_screen_toggle()  # Maximizes the figure window

# Plot efficient frontier
plt.plot(sig_p.transpose(), mu_p, '-k', linewidth=1.5)
plt.axhline(y=B/C, linestyle='--', color='k', linewidth=1.5)
plt.axvline(1/np.sqrt(C), linestyle='--', color='k', linewidth=1.5)

# Plot scatter points for stocks and MVP
plt.scatter(sig, z, color='k', marker='o', label='Stocks', alpha=1)
plt.scatter(1/np.sqrt(C), B/C, color='k', marker='o', label='MVP', alpha=1)

# Set title, limits, and labels
plt.title('Efficient Frontier', fontsize=16)
plt.xlim(0, 0.1)
plt.xlabel('Portfolio Risk', fontsize=16)
plt.ylabel('Portfolio Expected Return', fontsize=16)

# Add legend
plt.legend(["Efficient frontier", "B/C", "1/sqrt(C)", "Stocks", "MVP"], loc='upper left', fontsize=14)

# Save the figure
plt.savefig('MV.eps', format='eps')



# Calculate optimal weights
g = 1 / D * (A * (V1 @ np.ones((N, 1))) - B * (V1 @ z))
h = 1 / D * (C * (V1 @ z) - B * (V1 @ np.ones((N, 1))))

mup = np.linspace(0.001, 0.05, 10)  # Objective return of the portfolio
wp = g + h * mup

# Calculate global minimum variance portfolio weights
w_mvp = g + h * B / C

# Create a table to store weights
columns = ["GMVP"] + [str(round(mu * 100, 2)) + "%" for mu in mup]
index = Ret.columns[1:]
Weights = pd.DataFrame(index=index, columns=columns, dtype=float)

# Assign values to the table
Weights["GMVP"] = w_mvp.flatten()
Weights[columns[1:]] = wp

# Write weights to Excel file
Weights.to_excel("Weights_MV.xlsx")