clc 
clear

%% Data

%Import Stock Returns
[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
Factors=readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

Mkt=table2array(Factors(2:end,2))/100;
Rf=table2array(Factors(2:end,5))/100;

ExRet = table2array(Ret(:,2:end))-Rf;

%% Regress

%Rolling-Window Size
rw=40;

%Stocks
[nS, mS]=size(ExRet);
alpha=NaN(nS-rw+1,mS);
beta=NaN(nS-rw+1,mS);

for t=1:nS-rw+1 
    for i=1:mS
        mdl=fitlm(Mkt(t:t+rw-1,1),ExRet(t:t+rw-1,i));
        alpha(t,i)=table2array(mdl.Coefficients(1,1));
        beta(t,i)=table2array(mdl.Coefficients(2,1));
    end
end

%% Plots

Dates_for_plot = linspace(table2array(Ret(rw,1)),table2array(Ret(end,1)),20); % Sequence of 'n' 
Dates_for_plot=string(datetime(Dates_for_plot, "Format",'MMM-yyyy')); % Convert in date

%Stocks
%alpha
alpha_plot=figure(1);
plot(1:nS-rw+1,alpha(:,1),'Color','r','LineWidth',1.5); hold on
plot(1:nS-rw+1,alpha(:,2),'Color','k','LineWidth',1.5);
plot(1:nS-rw+1,alpha(:,3),'Color','b','LineWidth',1.5);
plot(1:nS-rw+1,alpha(:,4),'Color','g','LineWidth',1.5);
plot(1:nS-rw+1,alpha(:,5),'Color','m','LineWidth',1.5);
plot(1:nS-rw+1, zeros(nS-rw+1),'--k'); hold off
set(gca, 'Xtick', linspace(1,nS-rw+1,20),'FontSize',14)
xticklabels(Dates_for_plot)
xtickangle(45)
xlim([1 nS-rw+1])
legend(Ret(:,2:end).Properties.VariableNames)
title('Alpha')
saveas(alpha_plot, 'alpha_S','epsc')

%beta
beta_plot=figure(2);
plot(1:nS-rw+1,beta(:,1),'Color','r','LineWidth',1.5); hold on
plot(1:nS-rw+1,beta(:,2),'Color','k','LineWidth',1.5);
plot(1:nS-rw+1,beta(:,3),'Color','b','LineWidth',1.5);
plot(1:nS-rw+1,beta(:,4),'Color','g','LineWidth',1.5);
plot(1:nS-rw+1,beta(:,5),'Color','m','LineWidth',1.5);
plot(1:nS-rw+1, ones(nS-rw+1),'--k'); hold off
set(gca, 'Xtick', linspace(1,nS-rw+1,20),'FontSize',14)
xticklabels(Dates_for_plot)
xtickangle(45)
xlim([1 nS-rw+1])
legend(Ret(:,2:end).Properties.VariableNames)
title('Beta')
saveas(beta_plot, 'beta_S','epsc')

    
    