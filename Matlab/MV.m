% MV.m
% This script implements the Mean-Variance portfolio optimization for a set
% of assets and plots the efficient frontier.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clear
clc

% Read asset returns from an Excel file
[Ret] = readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');

% Extract return data (assuming it starts from the 2nd column)
R = table2array(Ret(:, 2:end));
N = size(R, 2); % Number of assets

% Calculate mean returns (expected returns)
z = (mean(R, 1))';
% Calculate standard deviations (risk)
sig = std(R, 1);
% Calculate the covariance matrix of returns
V = cov(R);
% Calculate the inverse of the covariance matrix
V1 = inv(V);

% Calculate the parameters A, B, C, and D for the efficient frontier
% A = z' * V1 * z
A = z' * V1 * z;
% B = 1' * V1 * z
B = z' * V1 * ones(N, 1);
% C = 1' * V1 * 1
C = ones(1, N) * V1 * ones(N, 1);
% D = A*C - B^2
D = A * C - B^2;

% Define a range of target portfolio expected returns for plotting the efficient frontier
mu_p = 0:0.0001:0.015;

% Calculate the variance of the efficient frontier for the range of expected returns
sig2_p = (1 / D) * (C * mu_p.^2 - 2 * B * mu_p + A);
% Calculate the standard deviation (risk) of the efficient frontier
sig_p = sqrt(sig2_p);

%% Plotting the Efficient Frontier

p = figure(1);
p.WindowState = 'maximized'; % Maximize the plot window

% Plot the efficient frontier (Risk vs. Expected Return)
plot(sig_p, mu_p, '-k', 'LineWidth', 1.5);
hold on

% Plot horizontal line at the expected return of the Global Minimum Variance Portfolio (GMVP)
line([0, 1 / sqrt(C)], [B / C, B / C], 'LineStyle', '--', 'Color', 'k', 'Linewidth', 1.5);
% Plot vertical line at the standard deviation of the Global Minimum Variance Portfolio (GMVP)
line([1 / sqrt(C), 1 / sqrt(C)], [0, B / C], 'LineStyle', '--', 'Color', 'k', 'Linewidth', 1.5);

% Scatter plot of individual stock risk and return
scatter(sig, z, "filled");
% Scatter plot of the Global Minimum Variance Portfolio (GMVP)
scatter(1 / sqrt(C), B / C, 'filled')

title('Efficient Frontier', 'FontSize', 16)
xlim([0, 0.1]) % Set the x-axis limits
xlabel('Portfolio Risk', 'FontSize', 16)
ylabel('Portfolio Expected Return', 'FontSize', 16)
legend("Efficient frontier", "B/C", "1/sqrt(C)", "Stocks", "MVP", 'Location', 'northwest', 'FontSize', 14)
saveas(p, "MV.eps", 'epsc') % Save the plot as an EPS file

%% Optimal weights Calculation

% Calculate parameters for portfolio weights formula
g = 1 / D * (A * (V1 * ones(N, 1)) - B * (V1) * z);
h = 1 / D * (C * (V1 * z) - B * (V1 * ones(N, 1)));

% Define a range of target expected returns for which to calculate optimal weights
mup = linspace(0.001, 0.05, 10); % This is the objective return of the portfolio (i.e. what we would like to obtain)
% Calculate optimal weights for the target expected returns
wp = g + h * mup;

%% Global minimum variance portfolio weights

% Calculate the weights for the Global Minimum Variance Portfolio (GMVP)
w_mvp = g + h * B / C;

% Create a table to store the calculated weights
Weights = table('Size', [5, 11], 'VariableTypes', repmat({'double'}, [1, 11]), 'VariableNames', ["GMVP", round(mup * 100, 2) + "%"], ...
    'RowNames', Ret.Properties.VariableNames(2:end));

% Populate the weights table with GMVP weights and optimal weights
Weights(:, :) = array2table([w_mvp, wp]);

% Write the weights table to an Excel file
writetable(Weights, "Weights_MV.xlsx", "FileType", "spreadsheet", "WriteVariableNames", true, ...
    "WriteRowNames", true);

% Plot the portfolio allocation (weights) for different target returns
weight_bar = figure(2);
bar(mup, wp, 'stacked')
legend(Ret.Properties.VariableNames(2:end), 'Location', 'northwest')
xticks(mup) % Set x-axis ticks to the target returns
xticklabels(round(mup * 100, 2) + "%") % Set x-axis labels to target returns in percentage
ylabel("Weight (%)")
xlabel("Expected return (%)")
title("Portfolio allocation")
weight_bar.Position = [100 100 800 400]; % Set figure position and size
saveas(weight_bar, "MV_p.eps", 'epsc') % Save the plot as an EPS file

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
