clc 
clear

%% Data

%Import Stock Returns
[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
Factors=readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');
Uncertainty=readtable('../Data/Uncertainty.xlsx', VariableNamingRule='preserve');

VIX=table2array(Uncertainty(2:end,2))./table2array(Uncertainty(1:end-1,2))-1;
EPU=table2array(Uncertainty(2:end,3))./table2array(Uncertainty(1:end-1,3))-1;

Mkt=table2array(Factors(2:end,2))/100;
Rf=table2array(Factors(:,5))/100;

ExRet = table2array(Ret(:,2:end))-Rf(2:end);

%% Testing the CAPM (One-Factor Model)

%Stocks
% VIX
[nS, mS]=size(ExRet);
alpha_VIX=NaN(2,mS);
beta1_VIX=NaN(2,mS);
beta2_VIX=NaN(2,mS);
R2_VIX=NaN(1,mS);


for i=1:mS
    mdl=fitlm([Mkt(:,1), VIX(:,1)],ExRet(:,i));
    alpha_VIX(1,i)=table2array(mdl.Coefficients(1,1));
    alpha_VIX(2,i)=table2array(mdl.Coefficients(1,3));
    beta1_VIX(1,i)=table2array(mdl.Coefficients(2,1));
    beta1_VIX(2,i)=table2array(mdl.Coefficients(2,3));
    beta2_VIX(1,i)=table2array(mdl.Coefficients(3,1));
    beta2_VIX(2,i)=table2array(mdl.Coefficients(3,3));
    R2_VIX(1,i)=mdl.Rsquared.Adjusted(1,1);
end


% EPU
alpha_EPU=NaN(2,mS);
beta1_EPU=NaN(2,mS);
beta2_EPU=NaN(2,mS);
R2_EPU=NaN(1,mS);

for i=1:mS
    mdl=fitlm([Mkt(:,1), EPU(:,1)],ExRet(:,i));
    alpha_EPU(1,i)=table2array(mdl.Coefficients(1,1));
    alpha_EPU(2,i)=table2array(mdl.Coefficients(1,3));
    beta1_EPU(1,i)=table2array(mdl.Coefficients(2,1));
    beta1_EPU(2,i)=table2array(mdl.Coefficients(2,3));
    beta2_EPU(1,i)=table2array(mdl.Coefficients(3,1));
    beta2_EPU(2,i)=table2array(mdl.Coefficients(3,3));
    R2_EPU(1,i)=mdl.Rsquared.Adjusted(1,1);
end

%% Display Results
Stocks_VIX=array2table(round([alpha_VIX;beta1_VIX;beta2_VIX;R2_VIX],5), 'VariableNames',Ret(:,2:end).Properties.VariableNames,...
    'RowNames',{'alpha','(alpha t-stat)', 'Mkt',' (Mkt t-stat)', 'VIX', '(VIX t-stat)','Adj. R2'})

Stocks_EPU=array2table(round([alpha_EPU;beta1_EPU;beta2_EPU;R2_EPU],5), 'VariableNames',Ret(:,2:end).Properties.VariableNames,...
    'RowNames',{'alpha','(alpha t-stat)', 'Mkt',' (Mkt t-stat)', 'EPU', '(EPU t-stat)','Adj. R2'})

writetable(Stocks_VIX,'CAPM_Stock_VIX.csv','WriteRowNames',true)
writetable(Stocks_EPU,'CAPM_Stock_EPU.csv','WriteRowNames',true)

