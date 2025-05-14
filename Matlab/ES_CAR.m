% ES_CAR.m
% This script performs an event study analysis to calculate Cumulative
% Abnormal Returns (CAR) and Cumulative Average Abnormal Returns (CAAR)
% around specific events, using a CAPM model for expected returns.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clc
clear

%% Settings
America = 1; % Flag: 1 to use only aviation disasters in America, 0 for all
cutoff = 150; % Cutoff for number of casualties to filter events
start_CAR = -1; % Start of the Cumulative Abnormal Return (CAR) window relative to event day (e.g., -1 means 1 day before event)
end_CAR = 5; % End of the CAR window relative to event day (e.g., 5 means 5 days after event)

start_CAPM = -250; % Start of the CAPM estimation window relative to event day
end_CAPM = -50; % End of the CAPM estimation window relative to event day

%% Import Data
% Import industry portfolio returns (assuming it includes Transport industry)
Portfolios = readtable("../Data/IndustryPortfolios.xlsx", VariableNamingRule="preserve");
% Import daily Fama-French Factors
Factors = readtable('../Data/FF_FactorsDaily.xlsx', VariableNamingRule='preserve');
% Import Events data
Events = readtable('../Data/Events.xlsx', VariableNamingRule='preserve');

% Filter events based on the America flag and casualty cutoff
if America == 0
    Events = Events(table2array(Events(:, 13)) > cutoff, 4);
else
    Events = Events(table2array(Events(:, 13)) > cutoff & table2array(Events(:, "Zone")) == 1, 4);
end

% Extract Market Excess Return (Mkt-RF) and Risk-Free Rate (RF) from Factors
% Convert from percentage to decimal
MktRet = table2array(Factors(1:end, 2)) / 100;
Rf = table2array(Factors(1:end, 5)) / 100;

% Extract Transport industry returns (assuming it's in the 14th column)
% Convert from percentage to decimal
Ret = table2array(Portfolios(1:end, 14)) ./ 100;
% Calculate Excess Returns for the Transport industry
ExRet = Ret - Rf;

% Convert event dates to datetime objects
Events = datetime(table2array(Events(:, 1)), 'InputFormat', 'dd-MM-yyyy');
% Convert returns dates to datetime objects
DatesReturn = datetime(table2array(Factors(:, 1)), 'ConvertFrom', 'yyyymmdd');

% Get the number of observations and stocks (assuming ExRet is the relevant data)
[n, nS] = size(ExRet); % n = number of observations, nS = number of stocks (should be 1 for Transport industry)

% Get the row index in the Returns dates that corresponds to each event date
[~, Dates] = ismember(Events, DatesReturn(:, 1));
nE = length(Events); % Number of events

% Handle cases where an event date might not be in the returns dates
for i=1:length(Events)
    while Dates(i)==0 % If Dates(i) is not a trading day, try next day
        Events(i,:)=Events(i,:)+caldays(1);
        [~,Dates(i)] = ismember(Events(i,:),DatesReturn(:,1));
    end
end


% Initialize matrices to store CAPM parameters (alpha and beta) for each event
alpha = NaN(nE, nS);
beta = NaN(nE, nS);

%% CAPM Estimation (First Stage)

% Loop through each event to estimate CAPM parameters in the estimation window
for i = 1:nE
    % Define the estimation window relative to the event date
    estimation_window_start = Dates(i, 1) + start_CAPM;
    estimation_window_end = Dates(i, 1) + end_CAPM;

    % Ensure the estimation window is within the bounds of the data
    if estimation_window_start >= 1 && estimation_window_end <= n
        % Fit CAPM model: ExRet = alpha + beta * MktRet + epsilon
        mdl = fitlm(MktRet(estimation_window_start:estimation_window_end, 1), ExRet(estimation_window_start:estimation_window_end, :));

        % Store the estimated alpha and beta for the current event
        alpha(i, :) = table2array(mdl.Coefficients(1, 1));
        beta(i, :) = table2array(mdl.Coefficients(2, 1));
    end
end

% Initialize matrix to store predicted returns
PredRet = NaN(abs(start_CAR) + abs(end_CAR) + 1, nE, nS);

% Calculate predicted returns for each stock and around each event, for the CAR window
for t = 1:abs(start_CAR) + abs(end_CAR) + 1
    for i = 1:nE
        for k = 1:nS
            % Check if the date for prediction is within the data bounds
            prediction_date_idx = Dates(i, 1) + start_CAR - 1 + t;
            if prediction_date_idx >= 1 && prediction_date_idx <= n
                % Predict return using estimated CAPM parameters and market return
                PredRet(t, i, k) = alpha(i, k) + beta(i, k) * (MktRet(prediction_date_idx, 1));
            end
        end
    end
end

% Preallocate matrix for observed returns around events
ObsRet_agg = NaN(abs(start_CAR) + abs(end_CAR) + 1, nE, nS);

% Get observed returns for each stock and around each event, for the CAR window
for i = 1:nE
    for t = 1:abs(start_CAR) + abs(end_CAR) + 1
        for k = 1:nS
            % Check if the date for observation is within the data bounds
            observation_date_idx = Dates(i, 1) + start_CAR - 1 + t;
            if observation_date_idx >= 1 && observation_date_idx <= n
                ObsRet_agg(t, i, k) = ExRet(observation_date_idx, k);
            end
        end
    end
end

% Get abnormal returns (Observed return - Predicted return)
AbnRet = ObsRet_agg - PredRet;

% Get cumulative abnormal returns (sum of abnormal returns over the CAR window)
CAR = cumsum(AbnRet, 1, 'omitnan'); % Cumsum along the time dimension, ignoring NaNs

% Get cumulative average abnormal returns (average of CARs across all events)
CAAR = squeeze(mean(CAR, 2, 'omitnan')) * 100; % Squeeze to remove singleton dimensions, mean across events, convert to percentage

%% Plotting CAAR

% Get dates vector for the plot (relative to event day 0)
date = (start_CAR:1:end_CAR)';

% Create a zero line for reference in the plot
zero = zeros(abs(start_CAR) + abs(end_CAR) + 1, 1);

% Plot the Cumulative Average Abnormal Returns (CAAR)
p = figure(1);
plot(date, CAAR, 'Color', 'b', 'LineWidth', 1.5); hold on
plot(date, zero, '--k'); hold off % Plot the zero line

% Configure plot appearance
title('Cumulative Average Abnormal Return (CAAR)', 'FontSize', 16)
xlabel('Days relative to event', 'FontSize', 16)
ylabel('CAAR (%)', 'FontSize', 16)
saveas(p, 'CAAR', 'epsc') % Save the plot as an EPS file

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895
