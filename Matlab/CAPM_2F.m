% CAPM_2F.m
% This script tests two-factor asset pricing models by adding either VIX
% changes or EPU changes as a second factor to the CAPM.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clc
clear

%% Data

% Import Stock Returns from an Excel file
[Ret] = readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
% Import Fama-French Factors from an Excel file
Factors = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');
% Import Uncertainty data (including VIX and EPU) from an Excel file
Uncertainty = readtable('../Data/Uncertainty.xlsx', VariableNamingRule='preserve');

% Calculate percentage change in VIX and EPU
% VIX change: (Current VIX / Previous VIX) - 1
VIX = table2array(Uncertainty(2:end, 2)) ./ table2array(Uncertainty(1:end - 1, 2)) - 1;
% EPU change: (Current EPU / Previous EPU) - 1
EPU = table2array(Uncertainty(2:end, 3)) ./ table2array(Uncertainty(1:end - 1, 3)) - 1;

% Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
% Convert from percentage to decimal, starting from the 2nd row of Factors
Mkt = table2array(Factors(2:end, 2)) / 100;
Rf = table2array(Factors(:, 5)) / 100;

% Calculate Excess Returns for stocks (Stock Return - Risk-Free Rate)
% Note: Using Rf from row 2 onwards to match the length of other series
ExRet = table2array(Ret(:, 2:end)) - Rf(2:end);

%% Testing the CAPM with VIX (Two-Factor Model)

% Get the size of the excess returns data
[nS, mS] = size(ExRet); % nS = number of observations, mS = number of stocks

% Initialize matrices to store regression results for CAPM + VIX model
alpha_VIX = NaN(2, mS); % Alpha coefficient and its t-statistic
beta1_VIX = NaN(2, mS); % Beta for Market factor and its t-statistic
beta2_VIX = NaN(2, mS); % Beta for VIX factor and its t-statistic
R2_VIX = NaN(1, mS); % Adjusted R-squared

% Loop through each stock to perform the regression with Market and VIX
for i = 1:mS
    % Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * VIX + epsilon_i
    mdl = fitlm([Mkt(:, 1), VIX(:, 1)], ExRet(:, i));

    % Extract and store the regression results
    alpha_VIX(1, i) = table2array(mdl.Coefficients(1, 1)); % Alpha
    alpha_VIX(2, i) = table2array(mdl.Coefficients(1, 3)); % Alpha t-stat
    beta1_VIX(1, i) = table2array(mdl.Coefficients(2, 1)); % Market Beta
    beta1_VIX(2, i) = table2array(mdl.Coefficients(2, 3)); % Market Beta t-stat
    beta2_VIX(1, i) = table2array(mdl.Coefficients(3, 1)); % VIX Beta
    beta2_VIX(2, i) = table2array(mdl.Coefficients(3, 3)); % VIX Beta t-stat
    R2_VIX(1, i) = mdl.Rsquared.Adjusted(1, 1); % Adjusted R-squared
end

%% Testing the CAPM with EPU (Two-Factor Model)

% Initialize matrices to store regression results for CAPM + EPU model
alpha_EPU = NaN(2, mS); % Alpha coefficient and its t-statistic
beta1_EPU = NaN(2, mS); % Beta for Market factor and its t-statistic
beta2_EPU = NaN(2, mS); % Beta for EPU factor and its t-statistic
R2_EPU = NaN(1, mS); % Adjusted R-squared

% Loop through each stock to perform the regression with Market and EPU
for i = 1:mS
    % Fit linear model: ExRet_i = alpha_i + beta1_i * Mkt + beta2_i * EPU + epsilon_i
    mdl = fitlm([Mkt(:, 1), EPU(:, 1)], ExRet(:, i));

    % Extract and store the regression results
    alpha_EPU(1, i) = table2array(mdl.Coefficients(1, 1)); % Alpha
    alpha_EPU(2, i) = table2array(mdl.Coefficients(1, 3)); % Alpha t-stat
    beta1_EPU(1, i) = table2array(mdl.Coefficients(2, 1)); % Market Beta
    beta1_EPU(2, i) = table2array(mdl.Coefficients(2, 3)); % Market Beta t-stat
    beta2_EPU(1, i) = table2array(mdl.Coefficients(3, 1)); % EPU Beta
    beta2_EPU(2, i) = table2array(mdl.Coefficients(3, 3)); % EPU Beta t-stat
    R2_EPU(1, i) = mdl.Rsquared.Adjusted(1, 1); % Adjusted R-squared
end

%% Display and Save Results

% Create a table for the CAPM + VIX results
Stocks_VIX = array2table(round([alpha_VIX; beta1_VIX; beta2_VIX; R2_VIX], 5), 'VariableNames', Ret(:, 2:end).Properties.VariableNames, ...
    'RowNames', {'alpha', '(alpha t-stat)', 'Mkt', ' (Mkt t-stat)', 'VIX', '(VIX t-stat)', 'Adj. R2'})

% Create a table for the CAPM + EPU results
Stocks_EPU = array2table(round([alpha_EPU; beta1_EPU; beta2_EPU; R2_EPU], 5), 'VariableNames', Ret(:, 2:end).Properties.VariableNames, ...
    'RowNames', {'alpha', '(alpha t-stat)', 'Mkt', ' (Mkt t-stat)', 'EPU', '(EPU t-stat)', 'Adj. R2'})

% Write the results tables to CSV files
writetable(Stocks_VIX, 'CAPM_Stock_VIX.csv', 'WriteRowNames', true)
writetable(Stocks_EPU, 'CAPM_Stock_EPU.csv', 'WriteRowNames', true)

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
