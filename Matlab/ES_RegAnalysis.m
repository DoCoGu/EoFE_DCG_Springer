% ES_RegAnalysis.m
% This script performs regression analysis on industry portfolio returns,
% including event studies and calendar effects (weekday and tax days).
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

clc
clear

% Setting parameters for the analysis
maxlag = 3; % Maximum lag for lagged returns in the regression
cutoff = 150; % Cutoff value for filtering events (based on column 13 in Events.xlsx)
fD = 3; % Number of future days to create event dummies for
America = 1; % Flag to filter events by Zone (1 for America, 0 for All)

%% Import Data

% Read industry portfolio returns from an Excel file
Ret = readtable("../Data/IndustryPortfolios.xlsx", VariableNamingRule="preserve");
% Extract dates from the first column and convert to datetime objects
DatesReturn = datetime(table2array(Ret(:, 1)), 'ConvertFrom', 'yyyyMMdd');

% Read events data from an Excel file
Events = readtable('../Data/Events.xlsx', VariableNamingRule="preserve");

% Filter events based on the cutoff and Zone (if America flag is set)
if America == 0
    % Filter for all events above the cutoff
    Events = Events(table2array(Events(:, 13)) > cutoff, 4);
else
    % Filter for events in Zone 1 (America) above the cutoff
    Events = Events(table2array(Events(:, 13)) > cutoff & table2array(Events(:, "Zone")) == 1, 4);
end

% Extract returns data (assuming it's in the 14th column)
Ret = table2array(Ret(:, 14));

% Convert event dates to datetime objects
Events = datetime(table2array(Events(:, 1)), "InputFormat", 'dd-MM-yyyy');

% Preallocate vectors for weekday dummies
dow = NaN(length(DatesReturn), 4);

% Create dummy variables for weekdays (Monday to Thursday)
for d = 1:4
    for i = 1:length(DatesReturn)
        if weekday(DatesReturn(i)) == d + 1 % Monday=2, Tuesday=3, etc.
            dow(i, d) = 1;
        else
            dow(i, d) = 0;
        end
    end
end

% Create events dummy variables
[Match] = ismember(DatesReturn, Events); % Find dates in Returns that match event dates

EDummy = zeros(length(DatesReturn), fD); % Initialize event dummy matrix

% Create dummies for event days and subsequent days (up to fD)
for i = 1:length(DatesReturn)
    if Match(i) == 1 % If the current date is an event date
        % Set dummies for the next fD days
        for fd_idx = 1:fD
            if (i + fd_idx) <= length(DatesReturn)
                EDummy(i + fd_idx, fd_idx) = 1;
            end
        end
    end
end
% Remove the first fD rows as they cannot have event dummies referring to previous events
EDummy(1:fD, :) = [];

% Create a dummy variable for tax days (Jan 1-5)
T = zeros(height(Ret), 1);
dday = day(DatesReturn);
dmonth = month(DatesReturn);

% Identify dates that are in January and between the 1st and 5th
T(dmonth == 1 & dday >= 1 & dday <= 5) = 1;

Res = []; % Initialize results matrix

% Create lagged returns matrix
RetLag = lagmatrix(Ret(:, 1), [1:maxlag]);

% Combine all independent variables into a single matrix X
X = [RetLag, dow, T, EDummy(:, 1:fD)];

% Fit a linear regression model
mdl = fitlm(X, Ret(:, 1));
coeff = round(table2array(mdl.Coefficients(:, 1)), 5); % Extract coefficients
se = table2array(mdl.Coefficients(:, 2)); % Extract standard errors

% Calculate t-statistics
Sig = coeff ./ se;
% Combine coefficients and t-statistics for display
res = [coeff'; Sig'];

% Define variable names for the results table
varnames = ["Mon", "Tue", "Wed", "Thu", "TaxDays"];
rString = "R_{t-" + [1:maxlag] + "}"; % Names for lagged returns
eString = "Event +" + [1:fD]; % Names for event dummies
var = ["Constant", rString, varnames, eString]; % All variable names

% Create a table to display the regression results
tab = table('Size', [2, size(coeff, 1)], 'VariableNames', var, 'VariableTypes', repmat({'double'}, [size(coeff, 1), 1]));
tab(:, :) = array2table(res); % Populate the table with results

% Write the results table to an Excel file based on the America flag
if America == 0
    writetable(tab, 'RegressionAnalysis_All.xlsx')
else
    writetable(tab, 'RegressionAnalysis_America.xlsx')
end

% This MATLAB script contains code examples and material from the book:
% "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
% You can find more information and download the book at: link.springer.com/book/9783031861895