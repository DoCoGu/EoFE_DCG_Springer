# MATLAB Scripts for "Essentials of Financial Economics"

This repository contains a collection of MATLAB scripts designed to accompany the book "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler. These scripts provide practical examples and implementations of key concepts discussed in the book, allowing readers to apply the theory and explore financial economic models.

## Scripts Overview

Here is a brief description of each MATLAB script included:

* **`ES_RegAnalysis.m`**: Performs regression analysis on industry portfolio returns, including event studies and calendar effects (weekday and tax days).
* **`getuncef.m`**: A helper function used to calculate the unconstrained efficient frontier for a given set of asset expected returns and covariance matrix.
* **`MV.m`**: Implements the Mean-Variance portfolio optimization for a set of assets and plots the efficient frontier.
* **`MV_rf.m`**: Implements Mean-Variance portfolio optimization with a risk-free asset, deriving and plotting the Capital Market Line (CML).
* **`BlackLitterman.m`**: Implements the Black-Litterman model for portfolio allocation, blending market equilibrium returns with investor views.
* **`CAPM.m`**: Tests the single-factor Capital Asset Pricing Model (CAPM) for a set of stocks by regressing excess stock returns on the market excess return.
* **`CAPM_2F.m`**: Tests two-factor asset pricing models by adding either VIX changes or EPU changes as a second factor to the CAPM.
* **`CAPM_3F.m`**: Tests the Fama-French three-factor model for a set of stocks by regressing excess stock returns on the market, SMB, and HML factors.
* **`CAPM_RW.m`**: Implements a rolling-window regression to estimate time-varying CAPM alpha and beta for a set of stocks.
* **`CS_EPU.m`**: Performs Fama-MacBeth cross-sectional regression analysis, including Economic Policy Uncertainty (EPU) and consumption growth as factors.
* **`CS_T.m`**: Performs Fama-MacBeth cross-sectional regression analysis, including temperature changes and consumption growth as factors.
* **`ES_CAR.m`**: Performs an event study analysis to calculate Cumulative Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR) around specific events, using a CAPM model for expected returns.

## Requirements

To run these scripts, you will need:

* MATLAB installed on your system.
* The necessary data files referenced in the scripts (e.g., `IndustryPortfolios.xlsx`, `FF_Factors.xlsx`, `Returns.xlsx`, etc.). The file paths in the scripts assume a specific directory structure (e.g., `../Data/`). You may need to adjust these paths based on where you store your data files.
* Depending on the specific script, you may need access to certain MATLAB toolboxes (e.g., Financial Toolbox, Statistics and Machine Learning Toolbox).

## Usage

1.  Ensure you have MATLAB installed and the required data files in the correct relative paths as specified in the scripts.
2.  Open the desired `.m` file in MATLAB.
3.  Run the script. The output may include plots, tables, or saved data files depending on the script's purpose.

## Book Reference

These MATLAB scripts contain code examples and material from the book:

"Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler