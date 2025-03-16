clc 
clear

%% Data

%Import Stock Returns
[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
Factors=readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

Mkt=table2array(Factors(2:end,2))/100;
Rf=table2array(Factors(2:end,5))/100;

ExRet = table2array(Ret(:,2:end))-Rf;

%% Testing the CAPM (One-Factor Model)

%Stocks
[nS, mS]=size(ExRet);
alpha_S=NaN(2,mS);
beta_S=NaN(2,mS);
R2_S=NaN(1,mS);
for i=1:mS
    mdl=fitlm(Mkt(:,1),ExRet(:,i));
    alpha_S(1,i)=table2array(mdl.Coefficients(1,1));
    alpha_S(2,i)=table2array(mdl.Coefficients(1,3));
    beta_S(1,i)=table2array(mdl.Coefficients(2,1));
    beta_S(2,i)=table2array(mdl.Coefficients(2,3));
    R2_S(1,i)=mdl.Rsquared.Adjusted(1,1);
end

%% Display Results
Results=array2table(round([alpha_S;beta_S;R2_S],5), 'VariableNames',Ret(:,2:end).Properties.VariableNames,...
    'RowNames',{'alpha','(alpha t-stat)', 'beta',' (beta t-stat)','Adj. R2'})


writetable(Results,'CAPM_Stock.csv','WriteRowNames',true)

