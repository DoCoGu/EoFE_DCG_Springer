% CS_EPU.m
% This script performs Fama-MacBeth cross-sectional regression analysis,
% including Economic Policy Uncertainty (EPU) and consumption growth as factors.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clc
clear

%% Data

% Import Portfolio Returns from an Excel file
[Ret] = readtable('../Data/Portfolios.xlsx', VariableNamingRule='preserve');
% Import Fama-French Factors from an Excel file
FFactors = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');
% Import Uncertainty data (including EPU) from an Excel file
Uncertainty = readtable('../Data/Uncertainty.xlsx', VariableNamingRule='preserve');
% Import Consumption data from an Excel file
Consumption = readtable('../Data/Consumption.xlsx', VariableNamingRule='preserve');

% Calculate percentage change in EPU and Consumption
% EPU change: (Current EPU / Previous EPU) - 1
EPU = table2array(Uncertainty(2:end, 3)) ./ table2array(Uncertainty(1:end - 1, 3)) - 1;
% Consumption change: (Current Consumption / Previous Consumption) - 1
C = table2array(Consumption(2:end, 2)) ./ table2array(Consumption(1:end - 1, 2)) - 1;

% Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
% Convert from percentage to decimal, starting from the 2nd row of Factors
Mkt = table2array(FFactors(2:end, 2)) / 100;
Rf = table2array(FFactors(2:end, 5)) / 100;

% Calculate Excess Returns for portfolios (Portfolio Return - Risk-Free Rate)
% Note: Using data from row 2 onwards for Returns and Rf to match lengths
ExRet = table2array(Ret(2:end, 2:end)) / 100 - Rf;

% Combine the factors into a single matrix
Factors = [Mkt, C, EPU]; % Factors: Market Excess Return, Consumption Change, EPU Change

%% First Stage Regression (Time-Series Regressions)

% Get the size of the excess returns and factors data
[n1, n2] = size(ExRet); % n1 = number of observations, n2 = number of portfolios
nF = size(Factors, 2); % nF = number of factors

% Initialize matrices to store results from the first stage
CoefAll = NaN(nF, n2); % Stores the factor betas for each portfolio
Res = NaN(n1, n2); % Stores the residuals from each regression

% Loop through each portfolio to perform time-series regression
for i = 1:n2
    % Regress portfolio excess returns on the factors
    % Model: ExRet_i = alpha_i + beta_Mkt*Mkt + beta_dC*dC + beta_EPU*EPU + epsilon_i
    [Coef, ~, Residual] = regress(ExRet(:, i), [ones(n1, 1), Factors]); % Add a column of ones for the intercept (alpha)

    % Store the factor betas (excluding the intercept)
    CoefAll(:, i) = Coef(2:end);
    % Store the residuals
    Res(:, i) = Residual;
end

% Calculate the variance-covariance matrix of the residuals
VarCovErr = cov(Res);

%% Second-Stage Regression (Cross-Sectional Regressions)

% Calculate the average excess return for each portfolio
MeanRet = mean(ExRet);
% Transpose the factor betas matrix for the cross-sectional regression
Betas = CoefAll';

% Perform cross-sectional regression of average returns on factor betas
% Model: MeanRet_i = lambda_0 + lambda_Mkt*beta_Mkt_i + lambda_dC*beta_dC_i + lambda_EPU*beta_EPU_i + eta_i
% Using 'Intercept', false assumes lambda_0 = 0
mdl = fitlm(Betas, MeanRet', 'Intercept', false);

% Extract and store the factor risk premia (Lambdas) and their standard errors
SE = table2array(mdl.Coefficients(:, 2)); % Standard Errors
Lambda = table2array(mdl.Coefficients(:, 1)); % Factor Risk Premia (Lambdas)
Tstat = Lambda ./ SE; % T-statistics (standard)

% Calculate Shanken-corrected standard errors
Sigma_f = cov(Factors); % Covariance matrix of factors
% Formula for Shanken-corrected covariance matrix of Lambdas
VarLam = ((Betas' * Betas)^-1 * Betas' * VarCovErr * Betas * ((Betas' * Betas)^-1) * ...
    (1 + Lambda' * ((Sigma_f)^-1) * Lambda) + Sigma_f) / n1;
SE_Shanken = sqrt(diag(VarLam)); % Shanken-corrected Standard Errors
Tstat_Shanken = Lambda ./ SE_Shanken; % Shanken-corrected T-statistics

% Perform T cross-sectional regressions (for Newey-West standard errors)
LambdaFull = NaN(n1, nF); % Matrix to store Lambdas from each cross-sectional regression
for j = 1:n1
    % Use returns from a single time period for the cross-sectional regression
    MeanRet = ExRet(j, :);
    % Fit cross-sectional model for the current time period
    mdl = fitlm(Betas, MeanRet', 'intercept', false);
    % Store the estimated Lambdas for this time period
    LambdaFull(j, :) = table2array(mdl.Coefficients(:, 1))';
end

% Calculate the mean of the Lambdas across all time periods
LambdaMean = mean(LambdaFull, 1);

% Calculate HAC-corrected (Newey-West) standard errors
X = ones(n1, 1); % Dummy predictor for HAC calculation
hac_cov = zeros(nF, 1); % Initialize vector for HAC variances
for k = 1:nF
    % Use the time series of estimated Lambdas for each factor
    y = LambdaFull(:, k);

    % Fit a dummy regression to use the HAC function
    mdl = fitlm(X, y, 'intercept', false);
    % Calculate HAC standard error for the current factor's Lambda time series
    hac_cov(k, 1) = hac(mdl, 'type', 'HAC', 'bandwidth', 2, 'weights', 'BT'); % Bandwidth and weights can be adjusted
end

SE_NW = sqrt(hac_cov); % Newey-West Standard Errors

Tstat_NW = LambdaMean' ./ SE_NW; % Newey-West T-statistics

%% Results

% Get the names of the portfolios
NamePort = Ret(:, 2:end).Properties.VariableNames;

% Create a table for the First Stage results (Factor Betas)
FirstStageReg = array2table(Betas, 'VariableNames', {'Mkt', 'dC', 'EPU'}, 'RowNames', NamePort)

% Create a table for the Second Stage results (Factor Risk Premia and T-stats)
SecondStage = array2table([Lambda'; Tstat'; Tstat_NW'; Tstat_Shanken'], 'VariableNames', ...
    {'Mkt', 'dC', 'EPU'}, 'RowNames', {'Lambda', 'tstat', 't-stat HAC', 't-stat Shanken'})

% Write the results tables to CSV files
writetable(FirstStageReg, 'FirstStage_EPU.csv', 'WriteRowNames', true)
writetable(SecondStage, 'SecondStage_EPU.csv', 'WriteRowNames', true)

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
