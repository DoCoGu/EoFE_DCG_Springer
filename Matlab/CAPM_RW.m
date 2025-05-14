% CAPM_RW.m
% This script implements a rolling-window regression to estimate time-varying
% CAPM alpha and beta for a set of stocks.
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

%% Rolling Regression

% Define the size of the rolling window
rw = 40; % Number of observations in each window

% Get the size of the excess returns data
[nS, mS] = size(ExRet); % nS = total number of observations, mS = number of stocks

% Initialize matrices to store rolling regression results
alpha = NaN(nS - rw + 1, mS); % Alpha coefficient for each window and stock
beta = NaN(nS - rw + 1, mS); % Beta coefficient for each window and stock

% Perform rolling window regression
% Outer loop iterates through the start of each window
for t = 1:nS - rw + 1
    % Inner loop iterates through each stock
    for i = 1:mS
        % Fit linear model for the current window: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
        % Use data from the current window (t to t+rw-1)
        mdl = fitlm(Mkt(t:t + rw - 1, 1), ExRet(t:t + rw - 1, i));

        % Extract and store the regression results for the current window
        alpha(t, i) = table2array(mdl.Coefficients(1, 1)); % Alpha coefficient
        beta(t, i) = table2array(mdl.Coefficients(2, 1)); % Beta coefficient
    end
end

%% Plots

% Prepare dates for plotting the rolling estimates
% Create a sequence of dates corresponding to the end of each rolling window
Dates_for_plot = linspace(table2array(Ret(rw, 1)), table2array(Ret(end, 1)), 20); % Select 20 points for x-axis ticks
Dates_for_plot = string(datetime(Dates_for_plot, "Format", 'MMM-yyyy')); % Convert to formatted date strings

% Plotting Rolling Alpha
alpha_plot = figure(1);
% Plot alpha for each stock over time
plot(1:nS - rw + 1, alpha(:, 1), 'Color', 'r', 'LineWidth', 1.5); hold on
plot(1:nS - rw + 1, alpha(:, 2), 'Color', 'k', 'LineWidth', 1.5);
plot(1:nS - rw + 1, alpha(:, 3), 'Color', 'b', 'LineWidth', 1.5);
plot(1:nS - rw + 1, alpha(:, 4), 'Color', 'g', 'LineWidth', 1.5);
plot(1:nS - rw + 1, alpha(:, 5), 'Color', 'm', 'LineWidth', 1.5);
% Add a horizontal line at zero for reference
plot(1:nS - rw + 1, zeros(nS - rw + 1), '--k');
hold off

% Configure plot appearance
set(gca, 'Xtick', linspace(1, nS - rw + 1, 20), 'FontSize', 14) % Set x-axis ticks
xticklabels(Dates_for_plot) % Set x-axis labels to dates
xtickangle(45) % Rotate x-axis labels for readability
xlim([1 nS - rw + 1]) % Set x-axis limits
legend(Ret(:, 2:end).Properties.VariableNames) % Add legend with stock names
title('Rolling Alpha (CAPM)')
saveas(alpha_plot, 'alpha_S', 'epsc') % Save the plot as an EPS file

% Plotting Rolling Beta
beta_plot = figure(2);
% Plot beta for each stock over time
plot(1:nS - rw + 1, beta(:, 1), 'Color', 'r', 'LineWidth', 1.5); hold on
plot(1:nS - rw + 1, beta(:, 2), 'Color', 'k', 'LineWidth', 1.5);
plot(1:nS - rw + 1, beta(:, 3), 'Color', 'b', 'LineWidth', 1.5);
plot(1:nS - rw + 1, beta(:, 4), 'Color', 'g', 'LineWidth', 1.5);
plot(1:nS - rw + 1, beta(:, 5), 'Color', 'm', 'LineWidth', 1.5);
% Add a horizontal line at one for reference (CAPM beta is expected to be 1 for the market)
plot(1:nS - rw + 1, ones(nS - rw + 1), '--k');
hold off

% Configure plot appearance
set(gca, 'Xtick', linspace(1, nS - rw + 1, 20), 'FontSize', 14) % Set x-axis ticks
xticklabels(Dates_for_plot) % Set x-axis labels to dates
xtickangle(45) % Rotate x-axis labels for readability
xlim([1 nS - rw + 1]) % Set x-axis limits
legend(Ret(:, 2:end).Properties.VariableNames) % Add legend with stock names
title('Rolling Beta (CAPM)')
saveas(beta_plot, 'beta_S', 'epsc') % Save the plot as an EPS file

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
