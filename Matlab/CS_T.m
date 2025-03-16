clc 
clear

%% Data

%Import Stock Returns
[Ret]=readtable('../Data/PortfoliosLong.xlsx', VariableNamingRule='preserve');
FFactors=readtable('../Data/FF_FactorsLong.xlsx', VariableNamingRule='preserve');
Temperature=readtable('../Data/Temperature.xlsx', VariableNamingRule='preserve');
Consumption=readtable('../Data/ConsumptionLong.xlsx', VariableNamingRule='preserve');


T=diff(log(table2array(Temperature(:,2))));

C=table2array(Consumption(2:end,2))./table2array(Consumption(1:end-1,2))-1;

Mkt=table2array(FFactors(2:end,2))/100;
Rf=table2array(FFactors(2:end,5))/100;

ExRet = table2array(Ret(2:end,2:end))/100-Rf;

Factors=[Mkt, C, T];

%% First Stage Regression
[n1,n2] = size(ExRet);
nF=size(Factors,2);

CoefAll=NaN(nF,n2);
Res=NaN(n1,n2);

for i=1:n2
    [Coef,~,Residual] = regress(ExRet(:,i),[ones(n1,1), Factors]);    
    CoefAll(:,i)=Coef(2:end);
    Res(:,i)=Residual;
end
VarCovErr=cov(Res);

%% Second-Stage Regression
MeanRet   = mean(ExRet);
Betas= CoefAll';

% Average returns
mdl= fitlm(Betas,MeanRet','Intercept',false);
SE = table2array(mdl.Coefficients(:,2));
Lambda = table2array(mdl.Coefficients(:,1));
Tstat = Lambda./SE;

% Shanken correction
Sigma_f      = cov(Factors);
VarLam     = ((Betas'*Betas)^-1 *Betas'*VarCovErr*Betas*((Betas'*Betas )^-1)*...
    (1 + Lambda'*((Sigma_f)^-1)*Lambda) + Sigma_f)/n1;
SE_Shanken         = sqrt(diag(VarLam));
Tstat_Shanken=Lambda./SE_Shanken;


% T cross-section regressions
LambdaFull=NaN(n1,nF);
for j = 1:n1
    MeanRet=ExRet(j,:);
    mdl= fitlm(Betas,MeanRet','intercept',false); 
    LambdaFull(j,:)=table2array(mdl.Coefficients(:,1))';
end

LambdaMean=mean(LambdaFull,1);

% HAC-corrected standard-errors
X = ones(n1, 1); 
hac_cov = zeros(nF, 1); 
for k = 1:nF
    y = LambdaFull(:, k); 
    
    mdl = fitlm(X, y,'intercept',false); 
    hac_cov(k, 1) = hac(mdl, 'type', 'HAC','bandwidth',2, 'weights','BT');
end

SE_NW = sqrt(hac_cov);

Tstat_NW=LambdaMean'./SE_NW;

%% Results
NamePort=Ret(:,2:end).Properties.VariableNames;
FirstStageReg=array2table(Betas,  'VariableNames',{'Mkt','dC','T'},'RowNames',NamePort)
SecondStage=array2table([Lambda';Tstat'; Tstat_NW'; Tstat_Shanken'],'VariableNames',...
    {'Mkt','dC','dT'},'RowNames',{'Lambda','t-stat', 't-stat HAC', 't-stat Shanken'})

writetable(FirstStageReg,'FirstStage_T.csv','WriteRowNames',true)
writetable(SecondStage,'SecondStage_T.csv','WriteRowNames',true)




