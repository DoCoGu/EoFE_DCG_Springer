# Scripts for "Essentials of Financial Economics"

This repository contains a comprehensive collection of scripts designed to complement the material presented in the book **"Essentials of Financial Economics"** by Michael Donadelli, Michele Costola, and Ivan Gufler. These scripts serve as practical tools, allowing readers and practitioners to implement the theoretical concepts and financial economic models discussed throughout the book using various popular programming languages.

The code provided here is directly linked to the examples and methodologies detailed in the textbook, offering a hands-on approach to understanding and applying financial economic principles. By running and modifying these scripts, users can gain deeper insights into topics such as asset pricing, portfolio theory, and event studies.

## Cross-Platform Availability

To ensure broad accessibility and cater to different preferences, the scripts are provided in four widely used programming languages for quantitative analysis:

* **MATLAB**
* **Python**
* **Julia**
* **R**

While the core logic and functionality are consistent across languages, there may be minor differences in implementation details or required libraries specific to each environment.

## Scripts Overview

The collection covers a range of fundamental topics in financial economics. The main scripts and their general scope include:

* **Asset Pricing Models:** Scripts for testing classic asset pricing models like the Capital Asset Pricing Model (CAPM) and the Fama-French three-factor model, including variations with additional factors like VIX or Economic Policy Uncertainty (EPU), and rolling window estimations.
* **Portfolio Theory:** Scripts for implementing Mean-Variance portfolio optimization, both with and without a risk-free asset, and the Black-Litterman model to incorporate investor views.
* **Event Studies:** Scripts designed to analyze the impact of specific events (such as aviation disasters) on asset returns, calculating abnormal returns and cumulative abnormal returns, and exploring event-driven investment strategies.
* **Regression Analysis:** Scripts for performing regression analysis to identify factors influencing asset returns, including calendar effects and event-specific impacts.

### Helper Function: `getuncef`

The `getuncef` (Get Unconstrained Efficient Frontier) script is a standalone function designed to calculate the characteristics (expected return, variance, and weights) of portfolios lying on the unconstrained efficient frontier. This function is not intended to be run directly as a main analysis script, but rather is called or included by other scripts, particularly those related to portfolio optimization (`MV`, `MV_rf`, `BlackLitterman`), to perform a specific calculation needed within their analysis.

## Script Mapping by Language

Below is a list of the scripts provided, categorized by programming language. The core functionality is similar for scripts with the same base name across different languages.

* **MATLAB (.m)**:
    * `CAPM.m`
    * `CAPM_2F.m`
    * `CAPM_3F.m`
    * `CAPM_RW.m`
    * `CS_EPU.m`
    * `CS_T.m`
    * `ES_CAR.m`
    * `ES_InvestmentStrategy.m`
    * `ES_RegAnalysis.m`
    * `MV.m`
    * `MV_rf.m`
    * `BlackLitterman.m`
    * `getuncef.m`

* **Python (.py)**:
    * `CAPM.py`
    * `CAPM_2F.py`
    * `CAPM_3F.py`
    * `CAPM_RW.py`
    * `CS_EPU.py`
    * `CS_T.py`
    * `ES_CAR.py`
    * `ES_InvestmentStrategy.py`
    * `ES_RegAnalysis.py`
    * `MV.py`
    * `MV_rf.py`
    * `BlackLitterman.py`
    * `getuncef.py`

* **Julia (.jl)**:
    * `CAPM.jl`
    * `CAPM_2F.jl`
    * `CAPM_3F.jl`
    * `CAPM_RW.jl`
    * `CS_EPU.jl`
    * `CS_T.jl`
    * `ES_CAR.jl`
    * `ES_InvestmentStrategy.jl`
    * `MV.jl`
    * `MV_rf.jl`
    * `BlackLitterman.jl`
    * `getuncef.jl`

* **R (.R)**:
    * `CAPM.R`
    * `CAPM_2F.R`
    * `CAPM_3F.R`
    * `CAPM_RW.R`
    * `CS_EPU.R`
    * `CS_T.R`
    * `ES_CAR.R`
    * `ES_InvestmentStrategy.R`
    * `ES_RegAnalysis.R`
    * `MV.R`
    * `MV_rf.R`
    * `BlackLitterman.R`
    * `getuncef.R`

## Data Structure

The scripts assume the presence of a `Data` folder located one directory level above the script folders (e.g., if your scripts are in `./Python`, the `Data` folder should be in `../Data`). This `Data` folder should contain the necessary Excel (`.xlsx`) data files referenced by the scripts (e.g., `Returns.xlsx`, `FF_Factors.xlsx`, `Events.xlsx`, etc.).

## Disclaimer and Contact

These scripts are provided as supplementary material to the book and are intended for educational and illustrative purposes. While efforts have been made to ensure their accuracy and consistency with the book's examples, they may contain bugs or errors.

The authors welcome feedback and bug reports. If you encounter any issues or have suggestions for improvement, please feel free to reach out to the authors:

* Michael Donadelli
* Michele Costola
* Ivan Gufler

*(Please refer to the book or publisher's website for contact information.)*

Using these scripts requires a basic understanding of the respective programming language and the concepts from the book. Users are encouraged to explore the code, experiment with different data, and adapt the scripts to their specific needs.

---