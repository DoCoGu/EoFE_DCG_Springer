% BlackLitterman.m
% This script implements the Black-Litterman model for portfolio allocation,
% blending market equilibrium returns with investor views.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clear
clc

%% Import and prepare data

% Read asset returns from an Excel file (assuming 5 assets)
ret = readtimetable('../Data/Returns.xlsx', VariableNamingRule='preserve');
% Read Fama-French factors from an Excel file
ff = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

% Add Market Excess Return and Risk-Free Rate to the returns timetable
ret.("Mkt-RF") = ff.("Mkt-RF")(2:end) / 100; % Market Excess Return (converted to decimal)
ret.rf = ff.RF(2:end) / 100; % Risk-Free Rate (converted to decimal)
clear mkt; % Clear potential variable named mkt

% Calculate excess returns for assets (asset return - risk-free rate)
exret = ret{:, 1:end - 2} - ret{:, end};
% Extract Market Excess Return
exmkt = ret{:, end - 1};

% Compute the equilibrium returns based on CAPM
[T, n] = size(exret); % T = number of observations, n = number of assets
beta = nan(2, n); % Initialize matrix to store CAPM betas

% Estimate CAPM beta for each asset by regressing excess asset returns on market excess return
for i = 1:n
    mdl = fitlm(exmkt, exret(:, i));
    beta(:, i) = table2array(mdl.Coefficients(2, 1)); % Store the beta coefficient
end

% Calculate the implied market equilibrium returns (Pi) using CAPM: Pi = beta * Market_Excess_Return
muExret = mean(exret); % Mean excess return of assets
Pi = transpose(beta(2, :) .* muExret); % Implied equilibrium excess returns; 

Sigma = cov(exret); % Covariance matrix of the excess returns
tau = 1 / T; % Uncertainty factor (Litterman and He, 1999) - often set to a small value

%% Define the Q and P matrices for investor views

% Q is the expected returns on the views, kx1 (k is the number of views)
Q = [0.04 / 12; ... % View 1: Expected return for a combination of assets
     0.02 / 12; ... % View 2: Expected return for another combination
     0.10 / 12]; % View 3: Expected return for a single asset

% P is the pick matrix, kxn (k views, n assets)
% Defines which assets are involved in each view and their weights
P = [1 0 -1 0 0; ... % View 1: Asset 1 - Asset 3
     0 -1 0 1 0; ... % View 2: Asset 4 - Asset 2
     0 0 0 0 1]; % View 3: Asset 5

% Omega is the covariance matrix of the views, kxk
% Represents the uncertainty in the investor's views
Omega = P * (tau * Sigma) * P';

%% Blend the Equilibrium returns Pi with the views

% Calculate the Black-Litterman covariance matrix (SigmaBL)
SigmaBL = inv(inv(tau * Sigma) + P' / Omega * P);
% Calculate the Black-Litterman expected returns (muBL)
muBL = SigmaBL * ((tau * Sigma) \ Pi + P' / (Omega) * Q);

% Create a table to compare mean excess returns, implied equilibrium returns (Pi), and Black-Litterman returns (muBL)
muTable = array2table(round([muExret; Pi'; muBL'], 4));
muTable.Properties.VariableNames = ret.Properties.VariableNames(1:end - 2);
muTable.Properties.RowNames = {'exret', 'Pi', 'muBL'};
writetable(muTable, 'muTable.csv', 'WriteRowNames', true); % Write the table to a CSV file

%% Implement the two efficient frontiers (Mean-Variance and Black-Litterman)

% Define a range of target expected returns for plotting the efficient frontiers
muR = 0:0.00001:0.015;
% Calculate the Mean-Variance efficient frontier using the mean excess returns
[~, OptSigma, ~, ~, w] = getuncef(muExret', Sigma, muR);
% Calculate the Black-Litterman efficient frontier using the blended returns and covariance
[~, OptSigmaBL, ~, ~, wBL] = getuncef(muBL, Sigma + SigmaBL, muR);

% Create a table to store the efficient frontier data (returns and standard deviations)
ef_table = array2table([muR; sqrt(OptSigma); sqrt(OptSigmaBL)]');
ef_table.Properties.VariableNames = {'r', 'sigma', 'sigmaBL'};

%% Plot the efficient frontiers

p = figure(1);
p.WindowState = 'maximized'; % Maximize the plot window

% Plot the Mean-Variance efficient frontier
plot(table2array(ef_table(:, 2)), table2array(ef_table(:, 1)), '-b', 'LineWidth', 1.5);
hold on
% Plot the Black-Litterman efficient frontier
plot(table2array(ef_table(:, 3)), table2array(ef_table(:, 1)), '-r', 'LineWidth', 1.5);
hold off

title('Efficient Frontier', 'FontSize', 16)
xlim([0, 0.1]) % Set x-axis limits
xlabel('Portfolio Risk', 'FontSize', 16)
ylabel('Portfolio Expected Return', 'FontSize', 16)
legend("M-V", "B-L", 'Location', 'northwest', 'FontSize', 14)
saveas(p, "BL.eps", 'epsc') % Save the plot as an EPS file

%% Optimal weights for a specific target return

% Create a table to store optimal weights for both models at a specific target return (e.g., 0.0060)
Weights = table('Size', [2, 5], 'VariableTypes', repmat({'double'}, [1, 5]), 'VariableNames', ret.Properties.VariableNames(1:5), ...
    'RowNames', ["M-V", "B-L"]);

% Find the index corresponding to the target return 0.0060 in the efficient frontier table
target_return_idx = table2array(ef_table(:, 1)) == 0.0060;
% Populate the weights table with optimal weights for the target return
Weights(:, :) = array2table([w(:, target_return_idx)'; ...
    wBL(:, target_return_idx)']);

% Write the weights table to an Excel file
writetable(Weights, "Weights_BL.xlsx", "FileType", "spreadsheet", "WriteVariableNames", true, ...
    "WriteRowNames", true);

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
