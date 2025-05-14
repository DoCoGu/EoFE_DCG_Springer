% getuncef.m
% This function calculates the unconstrained efficient frontier for a given
% set of asset expected returns and covariance matrix.
%
% Developed for: Essentials of Financial Economics
% Authors: Michael Donadelli, Michele Costola, Ivan Gufler
% Date: May 8, 2025
%%

function [muR, OptSigma, gmv, tradeoff, w] = getuncef(mu, S, muR)

    % Function: getuncef - get unconstrained efficient frontier
    
    % Get the number of assets
    k = size(S, 2);
    
    % Calculate the parameters A, B, C, and D from the portfolio theory
    % A = mu' * inv(S) * mu
    A = mu' * inv(S) * mu;
    % B = 1' * inv(S) * mu
    B = ones(k, 1)' * inv(S) * mu;
    % C = 1' * inv(S) * 1
    C = ones(k, 1)' * inv(S) * ones(k, 1);
    % D = A*C - B^2
    D = A * C - B^2;
    
    % Calculate the variance of the efficient frontier for given expected returns (muR)
    OptSigma = C .* muR.^2 - (2 * B) .* muR + A;
    OptSigma = OptSigma / D;
    
    % Calculate the expected return and variance of the Global Minimum Variance Portfolio (GMVP)
    gmv = [B / C; 1 / C]; % [Expected Return; Variance]
    
    % Calculate the tradeoff between risk and return (slope of the efficient frontier)
    tradeoff = [A / B; A / B^2]; % [Expected Return at tangency; Variance at tangency] - Note: This interpretation might vary slightly based on context
    
    % Calculate the weights of the efficient frontier portfolios for given expected returns (muR)
    % This uses the formula for portfolio weights on the efficient frontier
    d = A * inv(S) * ones(k, 1) - B * inv(S) * mu;
    E = C * inv(S) * mu - B * inv(S) * ones(k, 1);
    w = d / D + E / D .* muR;
    
    % This MATLAB script contains code examples and material from the book:
    % "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
    % You can find more information and download the book at: link.springer.com/book/9783031861895
    