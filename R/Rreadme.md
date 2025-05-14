# R Scripts for "Essentials of Financial Economics"

This repository contains a collection of R scripts designed to accompany the book "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler. These scripts provide practical examples and implementations of key concepts discussed in the book, allowing readers to apply the theory and explore financial economic models using R.

## Scripts Overview

Here is a brief description of each R script included:

* **`CAPM.R`**: Tests the single-factor Capital Asset Pricing Model (CAPM) for a set of assets by regressing excess asset returns on the market excess return.
* **`CAPM_2F.R`**: Tests two-factor asset pricing models by adding either VIX changes or EPU changes as a second factor to the CAPM.
* **`CAPM_3F.R`**: Tests the Fama-French three-factor model for a set of assets by regressing excess asset returns on the market, SMB, and HML factors.
* **`CAPM_RW.R`**: Implements a rolling-window regression to estimate time-varying CAPM alpha and beta for a set of assets and plots the results.
* **`CS_EPU.R`**: Performs Fama-MacBeth cross-sectional regression analysis, including Economic Policy Uncertainty (EPU) and consumption growth as factors, and calculates Newey-West and Shanken standard errors.
* **`CS_T.R`**: Performs Fama-MacBeth cross-sectional regression analysis, including temperature changes and consumption growth as factors, and calculates Newey-West and Shanken standard errors.
* **`ES_CAR.R`**: Performs an event study analysis to calculate Cumulative Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR) around specific events, using a CAPM model for expected returns.
* **`ES_InvestmentStrategy.R`**: Analyzes different investment strategies, including buy-and-hold and an event-driven strategy based on aviation disasters, and plots the cumulative returns.
* **`ES_RegAnalysis.R`**: Performs regression analysis on industry portfolio returns, including event studies and calendar effects (weekday and tax days).
* **`MV.R`**: Implements the Mean-Variance portfolio optimization for a set of risky assets and plots the efficient frontier, also calculating optimal portfolio weights.
* **`MV_rf.R`**: Implements Mean-Variance portfolio optimization with a risk-free asset, deriving and plotting the Capital Market Line (CML), and calculates optimal portfolio weights.
* **`BlackLitterman.R`**: Implements the Black-Litterman model for portfolio allocation, blending market equilibrium returns with investor views, and plots the resulting efficient frontiers and calculates optimal weights.
* **`getuncef.R`**: A helper function used to calculate the unconstrained efficient frontier for a given set of asset expected returns and covariance matrix. This function is used by other scripts like `MV.R`, `MV_rf.R`, and `BlackLitterman.R`.

## Requirements

To run these scripts, you will need:

* R installed on your system.
* The following R packages:
    * `readxl`
    * `dplyr`
    * `lubridate`
    * `ggplot2`
    * `sandwich`
    * `lmtest`
    * `stats`
    * `writexl`
    * `scales`
    * `PerformanceAnalytics` (Optional, but often useful)
* These packages can be installed in R using the `install.packages()` function. For example:
    ```r
    install.packages(c("readxl", "dplyr", "lubridate", "ggplot2", "sandwich", "lmtest", "stats", "writexl", "scales", "PerformanceAnalytics"))
    ```
* The necessary data files referenced in the scripts (e.g., `IndustryPortfolios.xlsx`, `FF_Factors.xlsx`, `Returns.xlsx`, etc.). The file paths in the scripts assume a specific directory structure (e.g., `../Data/`). You may need to adjust these paths based on where you store your data files.
* Ensure helper functions like `getuncef.R` are accessible by the scripts that source them (e.g., in the same directory or in R's search path).

## Usage

1.  Ensure you have R installed and the required packages and data files are in place.
2.  Open RStudio or your preferred R environment.
3.  Set your working directory to the location where you have saved the R scripts using `setwd("path/to/your/scripts")`.
4.  Run a script using the `source()` function, for example:
    ```r
    source("CAPM.R")
    ```
5.  The script may generate output in the console, create plots, or save results to files (e.g., CSV or Excel), depending on its functionality.

## Book Reference

These R scripts contain code examples and material from the book:

"Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.

You can find more information and download the book at: [link.springer.com/book/9783031861895](https://link.springer.com/book/9783031861895)

## Authors

* Michael Donadelli
* Michele Costola
* Ivan Gufler

---