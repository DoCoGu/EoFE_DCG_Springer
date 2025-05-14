% CAPM_3F.m
% This script tests the Fama-French three-factor model for a set of stocks
% by regressing excess stock returns on the market, SMB, and HML factors.
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
% Import Fama-French Factors (including Mkt-RF, SMB, HML, RF) from an Excel file
Factors = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

% Extract Fama-French factors (Mkt-RF, SMB, HML) and Risk-Free Rate (RF)
% Convert from percentage to decimal, starting from the 2nd row of Factors
FF = table2array(Factors(2:end, 2:4)) / 100; % Mkt-RF, SMB, HML
Rf = table2array(Factors(:, 5)) / 100; % RF

% Calculate Excess Returns for stocks (Stock Return - Risk-Free Rate)
% Note: Using Rf from row 2 onwards to match the length of other series
ExRet = table2array(Ret(:, 2:end)) - Rf(2:end);

%% Testing the CAPM (3-Factor Model)

% Get the size of the excess returns data
[nS, mS] = size(ExRet); % nS = number of observations, mS = number of stocks

% Initialize matrices to store regression results for the 3-Factor model
alpha = NaN(2, mS); % Alpha coefficient and its t-statistic
beta1 = NaN(2, mS); % Beta for Market factor and its t-statistic
beta2 = NaN(2, mS); % Beta for SMB factor and its t-statistic
beta3 = NaN(2, mS); % Beta for HML factor and its t-statistic
R2 = NaN(1, mS); % Adjusted R-squared

% Loop through each stock to perform the Fama-French 3-Factor regression
for i = 1:mS
    % Fit linear model: ExRet_i = alpha_i + beta1_i*Mkt + beta2_i*SMB + beta3_i*HML + epsilon_i
    mdl = fitlm(FF(:, 1:3), ExRet(:, i));

    % Extract and store the regression results
    alpha(1, i) = table2array(mdl.Coefficients(1, 1)); % Alpha
    alpha(2, i) = table2array(mdl.Coefficients(1, 3)); % Alpha t-stat
    beta1(1, i) = table2array(mdl.Coefficients(2, 1)); % Market Beta
    beta1(2, i) = table2array(mdl.Coefficients(2, 3)); % Market Beta t-stat
    beta2(1, i) = table2array(mdl.Coefficients(3, 1)); % SMB Beta
    beta2(2, i) = table2array(mdl.Coefficients(3, 3)); % SMB Beta t-stat
    beta3(1, i) = table2array(mdl.Coefficients(4, 1)); % HML Beta
    beta3(2, i) = table2array(mdl.Coefficients(4, 3)); % HML Beta t-stat
    R2(1, i) = mdl.Rsquared.Adjusted(1, 1); % Adjusted R-squared
end

%% Display and Save Results

% Create a table to display the regression results
Stocks = array2table(round([alpha; beta1; beta2; beta3; R2], 5), ...
    'VariableNames', Ret(:, 2:end).Properties.VariableNames, ...
    'RowNames', {'alpha', '(alpha t-stat)', 'beta_mkt', 'beta_mkt (t-stat)', ...
    'beta_smb', '(beta_smb t-stat)', 'beta_hml', '(beta_hml t-stat)', 'Adj. R2'})

% Write the results table to a CSV file
writetable(Stocks, 'CAPM_3F_Stock.csv', 'WriteRowNames', true)

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
