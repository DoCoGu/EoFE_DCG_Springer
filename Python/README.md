# Python Scripts for "Essentials of Financial Economics"

This repository contains a collection of Python scripts designed to accompany the book "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler. These scripts provide practical examples and implementations of key concepts discussed in the book, allowing readers to apply the theory and explore financial economic models using Python.

## Scripts Overview

Here is a brief description of each Python script included:

* **`CAPM.py`**: Tests the single-factor Capital Asset Pricing Model (CAPM) for a set of assets by regressing excess asset returns on the market excess return.
* **`CAPM_2F.py`**: Tests two-factor asset pricing models by adding either VIX changes or EPU changes as a second factor to the CAPM.
* **`CAPM_3F.py`**: Tests the Fama-French three-factor model for a set of assets by regressing excess asset returns on the market, SMB, and HML factors.
* **`CAPM_RW.py`**: Implements a rolling-window regression to estimate time-varying CAPM alpha and beta for a set of assets.
* **`BlackLitterman.py`**: Implements the Black-Litterman model for portfolio allocation, blending market equilibrium returns with investor views. (Requires `getuncef.py`)
* **`CS_EPU.py`**: Performs Fama-MacBeth cross-sectional regression analysis, including Economic Policy Uncertainty (EPU) and consumption growth as factors.
* **`CS_T.py`**: Performs Fama-MacBeth cross-sectional regression analysis, including temperature changes and consumption growth as factors.
* **`ES_CAR.py`**: Performs an event study analysis to calculate Cumulative Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR) around specific events, using a CAPM model for expected returns.
* **`ES_InvestmentStrategy.py`**: Analyzes different investment strategies, including buy-and-hold and an event-driven strategy based on aviation disasters.
* **`getuncef.py`**: A helper function used to calculate the unconstrained efficient frontier for a given set of asset expected returns and covariance matrix. This function is used by other scripts like `BlackLitterman.py`.

## Requirements

To run these scripts, you will need:

* Python installed on your system (Python 3.6 or later is recommended).
* The following Python libraries:
    * `pandas`
    * `numpy`
    * `matplotlib`
    * `statsmodels`
    * `openpyxl` (for reading `.xlsx` files)
* These libraries can typically be installed using pip:
    ```bash
    pip install pandas numpy matplotlib statsmodels openpyxl
    ```
* The necessary data files referenced in the scripts (e.g., `IndustryPortfolios.xlsx`, `FF_Factors.xlsx`, `Returns.xlsx`, etc.). The file paths in the scripts assume a specific directory structure (e.g., `../Data/`). You may need to adjust these paths based on where you store your data files.
* Ensure the `getuncef.py` file is accessible by the scripts that use it (e.g., in the same directory or in your Python path).

## Usage

1.  Ensure you have Python installed and the required libraries and data files are in place.
2.  Open your terminal or command prompt.
3.  Navigate to the directory where the Python scripts are saved.
4.  Run a script using the Python interpreter, for example:
    ```bash
    python CAPM.py
    ```
5.  The script may generate output in the console, create plots, or save results to files (e.g., CSV or Excel), depending on its functionality.

## Book Reference

These Python scripts contain code examples and material from the book:

"Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.

You can find more information and download the book at: [link.springer.com/book/9783031861895](https://link.springer.com/book/9783031861895)

## Authors

* Michael Donadelli
* Michele Costola
* Ivan Gufler

---
