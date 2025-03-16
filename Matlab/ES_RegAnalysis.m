clc
clear
% Setting
maxlag=3;
cutoff=150;
fD=3;
America=1;
%% Import Data

Ret = readtable("../Data/IndustryPortfolios.xlsx",VariableNamingRule="preserve");
DatesReturn=datetime(table2array(Ret(:,1)),'ConvertFrom','yyyyMMdd');

Events=readtable('../Data/Events.xlsx',VariableNamingRule="preserve");

if America==0
    Events=Events(table2array(Events(:,13))>cutoff,4);
else
    Events=Events(table2array(Events(:,13))>cutoff & table2array(Events(:,"Zone"))==1,4);
end



Ret = table2array(Ret(:,14));

% Convert dates in serial date 
Events=datetime(table2array(Events(:,1)),"InputFormat",'dd-MM-yyyy');


% Preallocate vectors of weekday dummies
dow=NaN(length(DatesReturn),4);

% Create dummy for weekdays
for d=1:4
for i=1:length(DatesReturn)
    if weekday(DatesReturn(i))==d+1
        dow(i,d)=1;
    else
        dow(i,d)=0;
    end
end
end


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

Res=[];
RetLag=lagmatrix(Ret(:,1),[1:maxlag]);

X=[RetLag,dow,T,EDummy(:,1:fD)];


mdl=fitlm(X,Ret(:,1));
coeff=round(table2array(mdl.Coefficients(:,1)),5);
se = table2array(mdl.Coefficients(:,2));


Sig=coeff./se;
res=[coeff';Sig'];


varnames=["Mon","Tue","Wed","Thu", "TaxDays"];
rString="R_{t-"+[1:maxlag]+"}";
eString="Event +"+[1:fD];
var=["Constant",rString,varnames,eString];

tab=table('size',[2,size(coeff,1)],'VariableNames',var, 'variabletypes', repmat({'double'},[size(coeff,1),1]));
tab(:,:)=array2table(res);

if America==0
    writetable(tab,'RegressionAnalysis_All.xlsx')
else
    writetable(tab,'RegressionAnalysis_America.xlsx')
end




