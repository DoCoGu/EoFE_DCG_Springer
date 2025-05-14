# Julia Scripts for "Essentials of Financial Economics"

This repository contains a collection of Julia scripts designed to accompany the book "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler. These scripts provide practical examples and implementations of key concepts discussed in the book, allowing readers to apply the theory and explore financial economic models using Julia.

## Scripts Overview

Here is a brief description of each Julia script included:

* **`CAPM.jl`**: Tests the single-factor Capital Asset Pricing Model (CAPM) for a set of assets by regressing excess asset returns on the market excess return.
* **`CAPM_2F.jl`**: Tests two-factor asset pricing models by adding either VIX changes or EPU changes as a second factor to the CAPM.
* **`CAPM_3F.jl`**: Tests the Fama-French three-factor model for a set of assets by regressing excess asset returns on the market, SMB, and HML factors.
* **`CAPM_RW.jl`**: Implements a rolling-window regression to estimate time-varying CAPM alpha and beta for a set of assets and plots the results.
* **`CS_EPU.jl`**: Performs Fama-MacBeth cross-sectional regression analysis, including Economic Policy Uncertainty (EPU) and consumption growth as factors, and calculates Newey-West and Shanken standard errors.
* **`CS_T.jl`**: Performs Fama-MacBeth cross-sectional regression analysis, including temperature changes and consumption growth as factors, and calculates Newey-West and Shanken standard errors.
* **`ES_CAR.jl`**: Performs an event study analysis to calculate Cumulative Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR) around specific events, using a CAPM model for expected returns.
* **`ES_InvestmentStrategy.jl`**: Analyzes different investment strategies, including buy-and-hold and an event-driven strategy based on aviation disasters, and plots the cumulative returns.
* **`MV.jl`**: Implements the Mean-Variance portfolio optimization for a set of risky assets and plots the efficient frontier, also calculating optimal portfolio weights.
* **`MV_rf.jl`**: Implements Mean-Variance portfolio optimization with a risk-free asset, deriving and plotting the Capital Market Line (CML), and calculates optimal portfolio weights.
* **`getuncef.jl`**: A helper function used to calculate the unconstrained efficient frontier for a given set of asset expected returns and covariance matrix. This function is used by other scripts like `BlackLitterman.jl` (though the `BlackLitterman.jl` script was not provided in Julia).

## Requirements

To run these scripts, you will need:

* Julia installed on your system.
* The following Julia packages:
    * `DataFrames`
    * `CSV`
    * `GLM`
    * `XLSX`
    * `Statistics`
    * `LinearAlgebra`
    * `Plots`
    * `Dates`
    * `StatsModels`
    * `ShiftedArrays`
    * `CovarianceMatrices`
* These packages can be installed using the Julia package manager. Open the Julia REPL and run:
    ```julia
    using Pkg
    Pkg.add(["DataFrames", "CSV", "GLM", "XLSX", "Statistics", "LinearAlgebra", "Plots", "Dates", "StatsModels", "ShiftedArrays", "CovarianceMatrices"])
    ```
* The necessary data files referenced in the scripts (e.g., `IndustryPortfolios.xlsx`, `FF_Factors.xlsx`, `Returns.xlsx`, etc.). The file paths in the scripts assume a specific directory structure (e.g., `../Data/`). You may need to adjust these paths based on where you store your data files.
* Ensure helper functions like `getuncef.jl` are accessible by the scripts that include them (e.g., in the same directory or in Julia's LOAD_PATH).

## Usage

1.  Ensure you have Julia installed and the required packages and data files are in place.
2.  Open the Julia REPL or your preferred Julia environment.
3.  Navigate to the directory where the Julia scripts are saved using `cd("path/to/your/scripts")`.
4.  Run a script using the `include()` function, for example:
    ```julia
    include("CAPM.jl")
    ```
5.  The script may generate output in the console, create plots, or save results to files (e.g., CSV or Excel), depending on its functionality.

## Book Reference

These Julia scripts contain code examples and material from the book:

"Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.

You can find more information and download the book at: [link.springer.com/book/9783031861895](https://link.springer.com/book/9783031861895)

## Authors

* Michael Donadelli
* Michele Costola
* Ivan Gufler

---