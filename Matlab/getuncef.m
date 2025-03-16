function [muR,OptSigma,gmv,tradeoff,w] = getuncef(mu,S,muR)

% ESGSF2021 - michele.costola@unive.it
% Matlab Code
% Function: getuncef - get uncontrained efficient frontier


k = size(S,2);

A = mu'*inv(S)*mu;
B = ones(k,1)'*inv(S)*mu;
C = ones(k,1)'*inv(S)*ones(k,1);
D = A*C-B^2;

OptSigma = C.*muR.^2-(2*B).*muR+A;
OptSigma = OptSigma/D;

gmv = [B/C;1/C];
tradeoff = [A/B;A/B^2];
d = A*inv(S)*ones(k,1)-B*inv(S)*mu;
E = C*inv(S)*mu-B*inv(S)*ones(k,1);
w = d/D+E/D.*muR;