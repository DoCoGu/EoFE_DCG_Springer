clc
clear 
% Setting
cutoff=150;
America=1;
%% Import Data

Ret = readtable("../Data/IndustryPortfolios.xlsx",VariableNamingRule="preserve");
Rf = readtable("../Data/FF_FactorsDaily.xlsx",VariableNamingRule="preserve");
Rf=Rf.RF;
DatesReturn=datetime(table2array(Ret(:,1)),'ConvertFrom','yyyyMMdd');

Events=readtable('../Data/Events.xlsx',VariableNamingRule="preserve");

if America==0
    Events=Events(table2array(Events(:,13))>cutoff,4);
else
    Events=Events(table2array(Events(:,13))>cutoff & table2array(Events(:,"Zone"))==1,4);
end



%Ret = table2array(Ret(2:end,2))./table2array(Ret(1:end-1,2))-1;
Ret = table2array(Ret(:,14));


% Convert dates in serial date 
Events=datetime(table2array(Events(:,1)),"InputFormat",'dd-MM-yyyy');


% Create events dummy
[Match]=ismember(DatesReturn,Events);


EDummy=zeros(length(DatesReturn),1);

for i=1:length(DatesReturn)
    if Match(i)==1
        EDummy(i+1,1)=1;
        EDummy(i+2,2)=1;
        EDummy(i+3,3)=1;
    else
        EDummy(i+1,1)=0;
        EDummy(i+2,2)=0;
        EDummy(i+3,3)=0;
    end
end
EDummy(1:3,:)=[];

T=zeros(height(Ret),1);
dday = day(DatesReturn);
dmonth = month(DatesReturn);

T(dmonth==1  & dday >= 1 & dday <= 5)=1;



%% Investment Strategy

% Buy and Hold strategy
BuyandHold =cumsum(Ret-Rf);

% Aviation disasters
AviationStrategy = Ret-Rf;
AviationStrategy(EDummy(:,1)==1)= ...
    Rf(EDummy(:,1)==1 )- ...
    Ret(EDummy(:,1)==1);

AviationStrategy=cumsum(AviationStrategy);

Inv_plot=figure(1);
plot(DatesReturn,BuyandHold,'-k',LineWidth=1.5);  hold on
plot(DatesReturn,AviationStrategy,'--b',LineWidth=1.5)
legend("Buy and Hold", "Aviation disasters Strategy","FontSize",16,Location="best")
ylabel('Cumulative return', 'FontSize', 18)
xlabel('Date', 'FontSize', 18)
set(gcf, 'Position', [0, 0, 5000, 5000]);
saveas(Inv_plot,'ES_InvStrategy.eps','epsc')
