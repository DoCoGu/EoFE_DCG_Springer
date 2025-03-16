clear
clc

[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');
Factors=readtable('../Data/FF_Factors.xlsx', VariableNamingRule='preserve');

Rf=table2array(Factors(2:end,5))/1200;

R=table2array(Ret(:,2:end));
N=size(R,2);

z = (mean(R,1))';
sig = std(R,1);
V = cov(R);
V1 = inv(V);

H= (z-ones(N,1)*mean(Rf))'*V1*(z-ones(N,1)*mean(Rf));
sqrtH=sqrt(H);

mu_p=linspace(mean(Rf),0.009,100);

A = z'*V1*z;
B = z'*V1*ones(N,1);
C = ones(1,N)*V1*ones(N,1);
D = A*C - B^2;

% Variance and standard deviation
sig2_p = 1/H*(mu_p-mean(Rf)).^2;
sig_p = sqrt(sig2_p);

mu_t=mean(Rf)+H/(ones(1,N)*V1*(z-ones(N,1)*mean(Rf)));
sig_t = sqrtH/(ones(1,N)*V1*(z-ones(N,1)*mean(Rf)));



A = z'*V1*z;
B = z'*V1*ones(N,1);
C = ones(1,N)*V1*ones(N,1);
D = A*C - B^2;
mu_p=linspace(0.001,0.008,100);

% Variance and standard deviation
sig2_p = (1/D)*(C*mu_p.^2 - 2*B*mu_p + A);
sig_pp = sqrt(sig2_p);
sig_pp(mu_p<0.005)=NaN;
x=[sig_pp;mu_p];

%% Plot
p=figure(2);
p.WindowState = 'maximized';
plot(sig_p, mean(Rf)+sqrtH.*sig_p,'-b','LineWidth',1.5); hold on
ax = gca;
ax.YAxis.Exponent=0;
%plot(sig_tp,mu_tp,'-k','LineWidth',1.5); hold on
plot(sig_pp,mu_p,'-k','LineWidth',1.5); hold on
scatter(sig_t, mu_t,"filled")
%scatter(sig,z, "filled");
title('Capital market line', 'FontSize', 16)
xlabel('Portfolio Risk', 'FontSize',16) 
ylabel('Portfolio Expected Return', 'FontSize',16) 
legend("Efficient frontier","Portfolio","Tangency port.","Autoupdate", "off",'FontSize',16,Location="northwest")
xlim([0,0.05])
line([0, sig_t], [mu_t, mu_t],'LineStyle','--','Color','k','Linewidth',1.5)
Ep_=mean(Rf)-sqrtH.*sig_p;
Ep_(Ep_<-0.002)=NaN;
plot(sig_p,Ep_ ,'-b','LineWidth',1.5);
line([sig_t, sig_t], [min(Ep_), mu_t],'LineStyle','--','Color','k','Linewidth',1.5)
%text(-0.003,mu_t,"$\sigma_T$",Interpreter="latex",FontSize=12)

saveas(p, "MV_rf.eps","epsc")
%% Optimal weights1
mup=linspace(0.001,0.05,10);
wp= V1*(z-ones(N,1)*mean(Rf))*(mup-mean(Rf))/H;


% Tangency portfolio weights
wT = (V1*(z-ones(N,1)*mean(Rf)))/(ones(N,1)'*V1*(z-ones(N,1)*mean(Rf)));

% Weight on the risk-free asset
rf= 1-sum([wT, wp]);

Weights=table('Size', [6,11],'VariableTypes',repmat({'double'},[1,11]),'VariableNames',["wT",round(mup*100,2)+"%"],...
    'RowNames',[Ret.Properties.VariableNames(2:end),"Rf"]);

Weights(:,:)=array2table(round([[wT,wp];rf],3));

writetable(Weights,"Weights_MVrf.xlsx","FileType","spreadsheet","WriteVariableNames",true,...
    "WriteRowNames",true);
