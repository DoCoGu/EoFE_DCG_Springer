clear
clc

[Ret]=readtable('../Data/Returns.xlsx', VariableNamingRule='preserve');

R=table2array(Ret(:,2:end));
N=size(R,2);

z = (mean(R,1))';
sig = std(R,1);
V = cov(R);
V1 = inv(V);

A = z'*V1*z;
B = z'*V1*ones(N,1);
C = ones(1,N)*V1*ones(N,1);
D = A*C - B^2;

mu_p=0:0.0001:0.015;

% Variance and standard deviation
sig2_p = (1/D)*(C*mu_p.^2 - 2*B*mu_p + A);
sig_p = sqrt(sig2_p);


%% Plot
p=figure(1);
p.WindowState = 'maximized';
plot(sig_p, mu_p,'-k','LineWidth',1.5); hold on
line([0,1/sqrt(C)], [B/C,B/C],'LineStyle','--','Color','k','Linewidth',1.5); 
line([1/sqrt(C), 1/sqrt(C)],[0, B/C],'LineStyle','--','Color','k','Linewidth',1.5); 
scatter(sig,z, "filled");
scatter(1/sqrt(C),B/C,'filled')
title('Efficient Frontier', 'FontSize', 16)
xlim([0,0.1])
xlabel('Portfolio Risk', 'FontSize', 16) 
ylabel('Portfolio Expected Return', 'FontSize', 16) 
legend("Efficient frontier", "B/C", "1/sqrt(C)", "Stocks", "MVP", 'Location','northwest','FontSize',14)
saveas(p, "MV.eps",'epsc')
%% Optimal weights

g= 1/D*(A*(V1*ones(N,1))-B*(V1)*z);
h=1/D*(C*(V1*z)-B*(V1*ones(N,1)));

mup = linspace(0.001, 0.05,10); % This is the objective return of the portfolio (i.e. what we would like to obtain)
wp=g+h*mup;



%% Global minimum variance portfolio weights

w_mvp = g + h *B/C;

Weights=table('Size', [5,11],'VariableTypes',repmat({'double'},[1,11]),'VariableNames',["GMVP", round(mup*100,2)+"%"],...
    'RowNames',Ret.Properties.VariableNames(2:end));

Weights(:,:)=array2table([w_mvp,wp]);


writetable(Weights,"Weights_MV.xlsx","FileType","spreadsheet","WriteVariableNames",true,...
    "WriteRowNames",true);

weight_bar = figure(2);
bar(mup, wp,'stacked')
legend(Ret.Properties.VariableNames(2:end), 'Location','northwest')
xticks(mup)
xticklabels(round(mup*100,2)+"%")
ylabel("Weight (%)")
xlabel("Expected return (%)")
title("Portfolio allocation")
weight_bar.Position = [100 100 800 400];
saveas(weight_bar, "MV_p.eps",'epsc')


