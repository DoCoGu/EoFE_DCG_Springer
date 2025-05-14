% MV_rf.m
% This script implements Mean-Variance portfolio optimization with a
% risk-free asset, deriving and plotting the Capital Market Line (CML).
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clear
clc

% Read asset returns from an Excel file
[Ret] = readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
% Read Fama-French factors from an Excel file
Factors = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

% Extract risk-free rate (assuming it's in the 5th column, starting from the 2nd row)
% Convert from annual percentage to monthly decimal
Rf = table2array(Factors(2:end, 5)) / 1200;

% Extract asset return data (assuming it starts from the 2nd column)
R = table2array(Ret(:, 2:end));
N = size(R, 2); % Number of assets

% Calculate mean returns
z = (mean(R, 1))';
% Calculate standard deviations
sig = std(R, 1);
% Calculate the covariance matrix of returns
V = cov(R);
% Calculate the inverse of the covariance matrix
V1 = inv(V);

% Calculate parameters for the efficient frontier with a risk-free asset
% H = (z - rf*1)' * inv(S) * (z - rf*1)
H = (z - ones(N, 1) * mean(Rf))' * V1 * (z - ones(N, 1) * mean(Rf));
% sqrt(H) is related to the slope of the Capital Market Line
sqrtH = sqrt(H);

% Define a range of target portfolio expected returns for plotting the CML
mu_p = linspace(mean(Rf), 0.009, 100);

% Calculate the variance of portfolios on the Capital Market Line
sig2_p = 1 / H * (mu_p - mean(Rf)).^2;
% Calculate the standard deviation (risk) of portfolios on the CML
sig_p = sqrt(sig2_p);

% Calculate the expected return and standard deviation of the Tangency Portfolio
mu_t = mean(Rf) + H / (ones(1, N) * V1 * (z - ones(N, 1) * mean(Rf)));
sig_t = sqrtH / (ones(1, N) * V1 * (z - ones(N, 1) * mean(Rf)));

% Recalculate parameters A, B, C, D for the risky asset efficient frontier (for plotting comparison)
A = z' * V1 * z;
B = z' * V1 * ones(N, 1);
C = ones(1, N) * V1 * ones(N, 1);
D = A * C - B^2;
% Define a range of target expected returns for the risky asset efficient frontier
mu_p = linspace(0.001, 0.008, 100);

% Calculate the variance of the risky asset efficient frontier
sig2_p = (1 / D) * (C * mu_p.^2 - 2 * B * mu_p + A);
sig_pp = sqrt(sig2_p);
% Set standard deviations below a certain threshold to NaN for cleaner plotting
sig_pp(mu_p < 0.005) = NaN;
x = [sig_pp; mu_p]; % Combine for potential use (though not directly used in plot)

%% Plotting the Capital Market Line and Efficient Frontier

p = figure(2);
p.WindowState = 'maximized'; % Maximize the plot window

% Plot the Capital Market Line (CML)
plot(sig_p, mean(Rf) + sqrtH .* sig_p, '-b', 'LineWidth', 1.5);
hold on
ax = gca;
ax.YAxis.Exponent = 0; % Display y-axis labels without scientific notation

% Plot the risky asset efficient frontier
plot(sig_pp, mu_p, '-k', 'LineWidth', 1.5);

% Scatter plot of the Tangency Portfolio
scatter(sig_t, mu_t, "filled")

title('Capital market line', 'FontSize', 16)
xlabel('Portfolio Risk', 'FontSize', 16)
ylabel('Portfolio Expected Return', 'FontSize', 16)
legend("Efficient frontier", "Portfolio", "Tangency port.", "Autoupdate", "off", 'FontSize', 16, Location="northwest")
xlim([0, 0.05]) % Set x-axis limits

% Add dashed lines to highlight the Tangency Portfolio's position
line([0, sig_t], [mu_t, mu_t], 'LineStyle', '--', 'Color', 'k', 'Linewidth', 1.5)
% Calculate the lower part of the CML (borrowing/lending at risk-free rate)
Ep_ = mean(Rf) - sqrtH .* sig_p;
Ep_(Ep_ < -0.002) = NaN; % Set values below a threshold to NaN for cleaner plotting
plot(sig_p, Ep_, '-b', 'LineWidth', 1.5);
line([sig_t, sig_t], [min(Ep_), mu_t], 'LineStyle', '--', 'Color', 'k', 'Linewidth', 1.5)
% text(-0.003,mu_t,"$\sigma_T$",Interpreter="latex",FontSize=12) % Example of adding text annotation

saveas(p, "MV_rf.eps", "epsc") % Save the plot as an EPS file

%% Optimal weights Calculation

% Define a range of target expected returns for which to calculate optimal weights
mup = linspace(0.001, 0.05, 10);
% Calculate optimal weights for portfolios on the CML (combination of risk-free and tangency)
wp = V1 * (z - ones(N, 1) * mean(Rf)) * (mup - mean(Rf)) / H;

% Calculate the weights for the Tangency Portfolio (risky assets only)
wT = (V1 * (z - ones(N, 1) * mean(Rf))) / (ones(N, 1)' * V1 * (z - ones(N, 1) * mean(Rf)));

% Calculate the weight allocated to the risk-free asset (1 - sum of risky asset weights)
rf = 1 - sum([wT, wp]);

% Create a table to store the calculated weights
Weights = table('Size', [6, 11], 'VariableTypes', repmat({'double'}, [1, 11]), 'VariableNames', ["wT", round(mup * 100, 2) + "%"], ...
    'RowNames', [Ret.Properties.VariableNames(2:end), "Rf"]);

% Populate the weights table with Tangency Portfolio weights and optimal CML portfolio weights (including risk-free)
Weights(:, :) = array2table(round([[wT, wp]; rf], 3));

% Write the weights table to an Excel file
writetable(Weights, "Weights_MVrf.xlsx", "FileType", "spreadsheet", "WriteVariableNames", true, ...
    "WriteRowNames", true);

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
