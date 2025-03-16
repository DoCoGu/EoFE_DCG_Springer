clc 
clear

%% Data

%Import Stock Returns
[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
Factors=readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

FF=table2array(Factors(2:end,2:4))/100;
Rf=table2array(Factors(:,5))/100;

ExRet = table2array(Ret(:,2:end))-Rf(2:end);

%% Testing the CAPM (3-Factor Model)

%Stocks
[nS, mS]=size(ExRet);
alpha=NaN(2,mS);
beta1=NaN(2,mS);
beta2=NaN(2,mS);
beta3=NaN(2,mS);
R2=NaN(1,mS);

for i=1:mS
    mdl=fitlm(FF(:,1:3),ExRet(:,i));
    alpha(1,i)=table2array(mdl.Coefficients(1,1));
    alpha(2,i)=table2array(mdl.Coefficients(1,3));
    beta1(1,i)=table2array(mdl.Coefficients(2,1));
    beta1(2,i)=table2array(mdl.Coefficients(2,3));
    beta2(1,i)=table2array(mdl.Coefficients(3,1));
    beta2(2,i)=table2array(mdl.Coefficients(3,3));
    beta3(1,i)=table2array(mdl.Coefficients(4,1));
    beta3(2,i)=table2array(mdl.Coefficients(4,3));
    R2(1,i)=mdl.Rsquared.Adjusted(1,1);
end

%% Display Results
Stocks=array2table(round([alpha;beta1;beta2;beta3;R2],5),...
    'VariableNames',Ret(:,2:end).Properties.VariableNames,...
    'RowNames',{'alpha','(alpha t-stat)', 'beta_mkt','beta_mkt (t-stat)',...
    'beta_smb','(beta_smb t-stat)', 'beta_hml','(beta_hml t-stat)','Adj. R2'})


writetable(Stocks,'CAPM_3F_Stock.csv','WriteRowNames',true)

