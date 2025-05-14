% CAPM.m
% This script tests the single-factor Capital Asset Pricing Model (CAPM)
% for a set of stocks by regressing excess stock returns on the market
% excess return.
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

% Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF)
% Convert from percentage to decimal, starting from the 2nd row of Factors
Mkt = table2array(Factors(2:end, 2)) / 100;
Rf = table2array(Factors(2:end, 5)) / 100;

% Calculate Excess Returns for stocks (Stock Return - Risk-Free Rate)
ExRet = table2array(Ret(:, 2:end)) - Rf;

%% Testing the CAPM (One-Factor Model)

% Get the size of the excess returns data
[nS, mS] = size(ExRet); % nS = number of observations, mS = number of stocks

% Initialize matrices to store regression results
alpha_S = NaN(2, mS); % Alpha coefficient and its t-statistic
beta_S = NaN(2, mS); % Beta coefficient and its t-statistic
R2_S = NaN(1, mS); % Adjusted R-squared

% Loop through each stock to perform the CAPM regression
for i = 1:mS
    % Fit linear model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
    mdl = fitlm(Mkt(:, 1), ExRet(:, i));

    % Extract and store the regression results
    alpha_S(1, i) = table2array(mdl.Coefficients(1, 1)); % Alpha coefficient
    alpha_S(2, i) = table2array(mdl.Coefficients(1, 3)); % Alpha t-statistic
    beta_S(1, i) = table2array(mdl.Coefficients(2, 1)); % Beta coefficient
    beta_S(2, i) = table2array(mdl.Coefficients(2, 3)); % Beta t-statistic
    R2_S(1, i) = mdl.Rsquared.Adjusted(1, 1); % Adjusted R-squared
end

%% Display and Save Results

% Combine all results into a single table for display and saving
Results = array2table(round([alpha_S; beta_S; R2_S], 5), 'VariableNames', Ret(:, 2:end).Properties.VariableNames, ...
    'RowNames', {'alpha', '(alpha t-stat)', 'beta', ' (beta t-stat)', 'Adj. R2'})

% Write the results table to a CSV file
writetable(Results, 'CAPM_Stock.csv', 'WriteRowNames', true)

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
