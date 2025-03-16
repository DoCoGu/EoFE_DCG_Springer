clc
clear
%% Settings
America=1; % Use only aviation disasters in America
cutoff=150; % Cutoff for number of causalties
start_CAR=-1;
end_CAR=5;

start_CAPM=-250;
end_CAPM=-50;

%% Import Data
Portfolios = readtable("../Data/IndustryPortfolios.xlsx",VariableNamingRule="preserve");
Factors=readtable('../Data/FF_FactorsDaily.xlsx',VariableNamingRule='preserve');
Events=readtable('../Data/Events.xlsx',VariableNamingRule='preserve');

if America==0
    Events=Events(table2array(Events(:,13))>cutoff,4);
else
    Events=Events(table2array(Events(:,13))>cutoff & table2array(Events(:,"Zone"))==1,4);
end


MktRet=table2array(Factors(1:end,2))/100;
Rf=table2array(Factors(1:end,5))/100;

Ret = table2array(Portfolios(1:end,14))./100; % Transport inudstry
ExRet = Ret- Rf;

% Convert dates in serial date 
Events=datetime(table2array(Events(:,1)),'InputFormat','dd-MM-yyyy');
DatesReturn=datetime(table2array(Factors(:,1)),'ConvertFrom','yyyymmdd');

[n,nS]=size(ExRet);

% Get event's date in returns dates
[~,Dates]=ismember(Events,DatesReturn(:,1));
for i=1:length(Events)
    while Dates(i)==0 % If Dates(i) is not a trading day, try next day
        Events(i,:)=Events(i,:)+caldays(1);
        [~,Dates(i)] = ismember(Events(i,:),DatesReturn(:,1));
    end
end


Dates=Dates(Dates>abs(start_CAPM));
Dates=unique(Dates);
nE=length(Dates);

%% Event Study
% Preallocate vectors for CAPM coefficients      
alpha=NaN(nE,nS);
beta=NaN(nE,nS);
% Compute the CAPM model for each stock and around each event
for i=1:nE
    for k=1:nS
        onefctmdl= fitlm(MktRet(Dates(i,1)+start_CAPM:Dates(i,1)+end_CAPM,1), ExRet(Dates(i,1)+start_CAPM:Dates(i,1)+end_CAPM,k));
        alpha(i,k)=table2array(onefctmdl.Coefficients(1,1));
        beta(i,k)=table2array(onefctmdl.Coefficients(2,1));
    end
end

  
    
% Preallocate vector of CAPM predicted returns
PredRet=NaN(abs(start_CAR)+abs(end_CAR)+1,nE,nS);
% Get predicted returns for each stocka and around each event, for the CARs
% period
for t=1:abs(start_CAR)+abs(end_CAR)+1
    for i=1:nE
        for k=1:nS
            PredRet(t,i,k)=alpha(i,k)+beta(i,k)*(MktRet(Dates(i,1)+start_CAR-1+t,1));
        end
    end
end


% Preallocate vector of observed returns 
ObsRet_agg=NaN(abs(start_CAR)+abs(end_CAR)+1,nE,nS);
% Get observed returns for each stocka and around each event, for the CARs
% period
for i=1:nE
    for t=1:abs(start_CAR)+abs(end_CAR)+1
        for k=1:nS
        ObsRet_agg(t,i,k)=ExRet(Dates(i,1)+start_CAR-1+t,k);
        end
    end
end
    
% Get abnormal returns (Observed ret. - Predicted ret.)
AbnRet=ObsRet_agg-PredRet;

% Get cumulative abnormal returns
CAR=cumsum(AbnRet,1,'omitnan');
   
CAAR=squeeze(mean(CAR,2,'omitnan'))*100;

%% Plot

% Get dates vector for plot
date=(start_CAR:1:end_CAR)';
for t=1:abs(start_CAR)+abs(end_CAR)+1
    date(t,1)=start_CAR-1+t;
end

zero=zeros(abs(start_CAR)+abs(end_CAR)+1,1);

CAAR_Plot=figure(1);
plot(date,CAAR,'-b',date, zero, '--k','LineWidth',2); hold on
line([0 0], [min(CAAR)-0.001 max(CAAR)+0.001],'Color','black','LineStyle','-','LineWidth', 2)
ylabel('CAR', 'FontSize', 18)
xlabel('Days Relative to Events', 'FontSize', 18)
xlim([start_CAR end_CAR])
ylim([min(CAAR)-0.001 max(CAAR)+0.001]);
legend('CAR')
a = get(gca,'XTickLabel');
set(gcf, 'Position', [0, 0, 1000, 1000]);
if America==0
    saveas(CAAR_Plot,'CAAR_All.eps','epsc')
else
    saveas(CAAR_Plot,'CAAR_America.eps','epsc')
end





