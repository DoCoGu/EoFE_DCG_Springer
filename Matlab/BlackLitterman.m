clear
clc

%% Import and prepare data
ret = readtimetable('../Data/Returns.xlsx', VariableNamingRule='preserve'); % 5 assets
ff = readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve'); % FamaFrench factors
ret.("Mkt-RF") = ff.("Mkt-RF")(2:end)/100;
ret.rf = ff.RF(2:end)/100;
clear mkt;

exret = ret{:,1:end-2}-ret{:,end}; 
exmkt = ret{:,end-1}; 
% Compute the equilibrium returns
[T,n] = size(exret); 
beta = nan(2,n);
for i = 1:n
    mdl = fitlm(exmkt,exret(:,i)); 
    beta(:,i) = table2array(mdl.Coefficients(2,1));
end

muExret = mean(exret);
Sigma = cov(exret); % Covariance matrix of the excess returns
Pi = transpose(beta(2,:).*muExret); % Market equilibrium returns
tau = 1/T; % Uncertainty factor (Litterman and He, 1999)

%% Define the Q and P
% Q is the expected returns on the views, nx1
Q = [0.04/12; ...
    0.02/12; ...
     0.10/12];

% P is the weighting matrix, kxn 
P = [1 0 -1 0 0; ... 
     0 -1 0 1 0; ... 
     0 0 0 0 1];

% Omega is the covariance matrix of the views, nxn
Omega = P*(tau*Sigma)*P';

% Blend the Equilibrium returns Pi with the views
SigmaBL = inv(inv(tau*Sigma)+P'/Omega*P);
muBL = SigmaBL*((tau*Sigma)\Pi+P'/(Omega)*Q);

% Table of muExret, Pi, muBL
muTable = array2table(round([muExret;Pi';muBL'],4));
muTable.Properties.VariableNames = ret.Properties.VariableNames(1:end-2);
muTable.Properties.RowNames = {'exret','Pi','muBL'};
writetable(muTable,'muTable.csv','WriteRowNames',true);

% Implement the two efficient frontiers
muR = 0:0.00001:0.015;
[~,OptSigma,~,~,w]  = getuncef(muExret',Sigma,muR); 
[~,OptSigmaBL,~,~,wBL] = getuncef(muBL,Sigma+SigmaBL,muR); 

ef_table = array2table([muR;sqrt(OptSigma);sqrt(OptSigmaBL)]');
ef_table.Properties.VariableNames = {'r','sigma','sigmaBL'};

%% Plot the efficient frontiers
p=figure(1);
p.WindowState = 'maximized';
plot(table2array(ef_table(:,2)), table2array(ef_table(:,1)),'-b','LineWidth',1.5); hold on
plot(table2array(ef_table(:,3)), table2array(ef_table(:,1)),'-r','LineWidth',1.5); hold on
title('Efficient Frontier', 'FontSize', 16)
xlim([0,0.1])
xlabel('Portfolio Risk', 'FontSize', 16) 
ylabel('Portfolio Expected Return', 'FontSize', 16) 
legend("M-V", "B-L", 'Location','northwest','FontSize',14)
saveas(p, "BL.eps",'epsc')

%% Optimal weights

Weights=table('Size', [2,5],'VariableTypes',repmat({'double'},[1,5]),'VariableNames',ret.Properties.VariableNames(1:5),...
    'RowNames',["M-V", "B-L"]);

Weights(:,:) = array2table([w(:,table2array(ef_table(:,1)==0.0060))'; ...
    wBL(:,table2array(ef_table(:,1)==0.0060))']);

writetable(Weights,"Weights_BL.xlsx","FileType","spreadsheet","WriteVariableNames",true,...
    "WriteRowNames",true);
