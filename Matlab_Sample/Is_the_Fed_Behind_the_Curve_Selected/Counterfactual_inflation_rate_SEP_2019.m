tic;

clear all
close all
warning off

addpath lib

starting_year=76; % Choose a starting year for iteration, 76 = 2019Q1

lambda=1; % Choose lambda

TsimEY=200; % Number of draws for bands sim
TsimX=200;

dum_ZLBcstr=1; % Impose ZLB constraint

% Max horizon of forecast
H=16; % quarters
x=(1:H)'; %time for plots

rho=0.9; % Set var-cov matrix for simulations for EY

printfig=0; % =1 to store figures
saveres=0; % =1 to store results
   
name{1}='unemployment'; % Import forecasts
name{2}='pi';

%% Load IRs estimated from BVAR 

if dum_ZLBcstr==1
    load EstimationIR_MYrep/IRs_BVAR_FF4_MYrep.mat
else
    load EstimationIR_MYrep/IRs_BVAR_FF4_ShadowFFR.mat
end

IR_FF4=10*mIRF_eps; % Number: dim(mIRF_eps)=(Ndraws,Nvar,Hhor) with vars= MP shock, ur, pi, ffr

IR_FF4=[IR_FF4(:,[1 2 4],:)]; % Re-organize IRs first so they all have same ordering


load EstimationIR_MYrep/IRs_BVAR_OR10_MYrep.mat mIRF_eps

IR_OR10=10*mIRF_eps; % NB: dim(mIRF_eps)=(Ndraws,Nvar,Hhor) with vars= MP shock, ur, pi, ffr

IR_OR10=[IR_OR10(:,[1 2 4],:)]; % Re-organize IRs first so they all have same ordering

IR=[IR_FF4 IR_OR10]; % Put all IRs togeteher

%% GET Raw EY0

for i=1:length(name)
    data=xlsread('fomc_projections_MYrep.xlsx',name{i},'A2:AK122');
    
    colNames ={'year','month','yearn','meeting_date','Actual0','Actual0_RT','median_t0','median_t','median_t1','median_t2','median_t3','median_lr','RealTime_LR','cent_low_t','cent_high_t','cent_low_t1','cent_high_t1','cent_low_t2','cent_high_t2','cent_low_t3','cent_high_t3','cent_low_lr','cent_high_lr','range_low_t','range_high_t','range_low_t1','range_high_t1','range_low_t2','range_high_t2','range_low_t3','range_high_t3','range_low_lr','range_high_lr','ffr','nber','T10YFF','shffr'};
    sTable = array2table(data,'VariableNames',colNames);
    
    % Get rid of obs with only one forecast over 81-07, so only keep 81-07 forecasts as of July. Keep everything post 07
    ind=find(sTable.year<2007 & sTable.month==1);
    sTable(ind,:)=[];
    
    year=sTable.year;
    month=sTable.month;
    yearn=year+(month-1)/12;
    
    quarter=[floor(((month-1)/3))+1];
    for iQ=1:length(quarter)
        if quarter(iQ)==1, Qinit(iQ)=5;end;
        if quarter(iQ)==2, Qinit(iQ)=4;end;
        if quarter(iQ)==3, Qinit(iQ)=3;end;
        if quarter(iQ)==4, Qinit(iQ)=2;end;
    end
    
    ffr=sTable.ffr;

    if i == 1
        sTable_UR = sTable;
    else
        sTable_PI = sTable;
    end
    
    % Fill up SEP's empty value with real-time longrun, as median_lr
    RTlr=sTable.median_lr;
    RTlr(isnan(sTable.median_lr))=sTable.RealTime_LR(isnan(sTable.median_lr));
    
    % Forecast data are as follows: last Q4 (backcast data), this year's Q4, next year's Q4, t+2yr, t+3yr and long run (t+5yr say)
    dataEu=[ sTable.median_t sTable.median_t1 sTable.median_t2 sTable.median_t3 RTlr];

  
    
    % Load sigma computed from normal distribution fitted on sep forecast range and central tendency
    load SIGMA_SEP_MYrep.mat SIGMA
    dataEuSig=SIGMA(:,:,i);
    
    %------------------------Revise backcast data-------------------------%
    
    backdata  =  xlsread('Counterfactual_Objectives.xlsx', 'A2:H178');
    colNames = {'yearn','year','quarter','UR','PI','FFR'};
    bTable = array2table(backdata,'VariableNames',colNames);
    bTable_raw = bTable;

    if i == 1
        target_variable = 'UR';
    else
        target_variable = 'PI';
    end
    
    sTable.backcast = NaN(height(sTable), 1); % Create an empty matrix for backcast data
    
    % Fill the line with last year Q4's data
    for j = 1:height(sTable)
        current_year = sTable.year(j);
        idx = find(bTable.year == current_year - 1 & bTable.quarter == 4);
        if ~isempty(idx)
            sTable.backcast(j) = bTable.(target_variable)(idx);
        end
    end

    dataEub=sTable.backcast; 
    dataEubActual=sTable.backcast;
    %---------------------------------------------------------------------%

    Eu=NaN(length(dataEu),H); 
    EuSig=NaN(length(dataEu),H);
    
    % Start filling the forecast at different horizon depending on month of FOMC (using "j"), and similarly to fill the backcast ("jb")
    for kk=1:length(dataEu)
        if month(kk)==1, j=3; jb=0; end;
        if month(kk)==3 | month(kk)==4, j=3; jb=0; end;
        if month(kk)==6, j=2; jb=1; end;
        if month(kk)==9, j=1; jb=2; end;
        if month(kk)==12 | month(kk)==11, j=0; jb=3; end;
        
        for ii=1:4
            Eu(kk,j+4*(ii-1)+1)=dataEu(kk,ii);
            EuSig(kk,j+4*(ii-1)+1)=dataEuSig(kk,ii);
        end
        
        %fill the long-run value for year 5
        Eu(kk,end)=dataEu(kk,end);
        EuSig(kk,end)=dataEuSig(kk,end);
        
        % add backcast to help with interpolation: 
        % use actual backcast if real time value is not available:
        if isnan(dataEub(kk))~=1
            ExtEu{kk}=[dataEub(kk) NaN(1,jb) Eu(kk,:)];
        else
            ExtEu{kk}=[dataEubActual(kk) NaN(1,jb) Eu(kk,:)];
        end
        ExtEuSig{kk}=[0 NaN(1,jb) EuSig(kk,:)];
        
        %interpolate missing nans
        temp1=naninterp(ExtEu{kk});
        temp2=naninterp(ExtEuSig{kk});
        
        %remove backcast to start at t=0
        Eui(kk,:)=temp1(2+jb:end);
        EuiSig(kk,:)=temp2(2+jb:end);
        
        
        %Forecast for the mandates in deviations from their targets (the real-time estimates of the long-run values)
        Eui_gap(kk,:)=Eui(kk,:)-dataEu(kk,end);
        
        %date and horizon associated with forecast
        SEP_hor(:,kk)=(year(kk)+((quarter(kk)-1)/4):1/4:year(kk)+((quarter(kk)-1)/4)+(H-1)/4)';
        
    end
    
    %store output
    Ey(:,:,i)=Eui;
    Ey_gap(:,:,i)=Eui_gap;
    EySig(:,:,i)=EuiSig;
    
end

%fill the nans with zeros
Ey(isnan(Ey))=0; Ey_gap(isnan(Ey_gap))=0;

%------------------------
%Import ffr path forecast

data_tp=xlsread('fomc_projections_MYrep.xlsx','ffr','A2:AJ122');

colNames_tp ={'year','month','yearn','meeting_date','Actual','RealTime','median_t0','median_t','median_t1','median_t2','median_t3','median_lr','RealTime_LR','cent_low_t','cent_high_t','cent_low_t1','cent_high_t1','cent_low_t2','cent_high_t2','cent_low_t3','cent_high_t3','cent_low_lr','cent_high_lr','range_low_t','range_high_t','range_low_t1','range_high_t1','range_low_t2','range_high_t2','range_low_t3','range_high_t3','range_low_lr','range_high_lr'};
sTable_tp = array2table(data_tp,'VariableNames',colNames_tp);

%Get rid of obs with only one forecast over 81-07, so only keep 81-07 forecasts as of July. Keep everything post 07
ind=find(sTable_tp.year<2007 & sTable_tp.month==1);
sTable_tp(ind,:)=[];


%As long-run value, use SEP and if not available, use real-time long-run estimate
RTlr=sTable_tp.median_lr;
RTlr(isnan(sTable_tp.median_lr))=sTable_tp.RealTime_LR(isnan(sTable_tp.median_lr));

%forecast data are as follows: last Q4, this year's Q4, next year's Q4, t+2yr, t+3yr and long run (t+5yr say)
dataEu=[ sTable_tp.median_t sTable_tp.median_t1 sTable_tp.median_t2 sTable_tp.median_t3 RTlr];

%backcast data

%-------------------------------Revise FFR forecast-----------------------%
backdata  =  xlsread('Counterfactual_Objectives.xlsx', 'A2:H178');
colNames = {'yearn','year','quarter','UR','PI','FFR'};
bTable = array2table(backdata,'VariableNames',colNames);


sTable_tp.backcast = NaN(height(sTable_tp), 1);

% Fill the line with last year Q4's data
for j = 1:height(sTable_tp)
    current_year = sTable_tp.year(j);
    idx = find(bTable.year == current_year - 1 & bTable.quarter == 4);
    if ~isempty(idx)
        sTable_tp.backcast(j) = bTable.FFR(idx);
    end
end

% Fine the indexes before 2015M6, and set them as NaN
idx = sTable_tp.year < 2015 | (sTable_tp.year == 2015 & sTable_tp.month <= 6);
sTable_tp.backcast(idx) = NaN;

dataEub=[]; 
dataEubActual=sTable_tp.backcast;
%-------------------------------------------------------------------------%

dataEU_now=sTable_tp.median_t0;

Eu=NaN(length(dataEu),H);
Eu_lb=NaN(length(dataEu),H);Eu_ub=NaN(length(dataEu),H);
Eu_lb2=NaN(length(dataEu),H);Eu_ub2=NaN(length(dataEu),H);
EuSig=NaN(length(dataEu),H);

year_tp=sTable_tp.year;
yearn_tp=sTable_tp.yearn;
month_tp=sTable_tp.month;

Eui=[];

%start filling the forecast at different horizon depending on month of FOMC (using "j"), and similarly to fill the backcast ("jb")
for kk=1:length(dataEu)
    if month_tp(kk)==1, j=3; jb=0; end;
    if month_tp(kk)==3 | month_tp(kk)==4, j=3; jb=0; end;
    if month_tp(kk)==6, j=2; jb=1; end;
    if month_tp(kk)==9, j=1; jb=2; end;
    if month_tp(kk)==12 | month_tp(kk)==11, j=0; jb=3; end;
    
    for ii=1:4
        Eu(kk,j+4*(ii-1)+1)=dataEu(kk,ii);
     end
    
    %fill the long-run value for year 5
    Eu(kk,end)=dataEu(kk,end);
    
    Eui(kk,:)=Eu(kk,:);
    
    %now add current reading (since it could be different from linear interpolation)
    Eui(kk,1)=dataEU_now(kk);
    
    Eui(kk,:)=naninterp( Eui(kk,:));
    
end

%store output
Ey(:,:,3)=Eui;
Ey_gap(:,:,3)=Eui;
EySig(:,:,3)=NaN(length(dataEu),H);

Ey_gap_raw = Ey_gap;

%%——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————%%


%% Calculate OPP and its confidence bands

%permute IR dims so it's easier to manipulate:
%Hhor,Nvar,Ndraw
IR=permute(IR,[3 2 1]);

%Only keep the first H horizons
Ey=Ey(:,1:H,:); Ey_gap=Ey_gap(:,1:H,:); EySig=EySig(:,1:H,:);
IR=IR(1:H,:,:);

%draw from IR
INDrnd=randsample(1:size(IR,3),TsimX,'true');

Xsim=IR(:,:,INDrnd);

X1sim=squeeze([IR(:,1,INDrnd); IR(:,2,INDrnd)]); %put ir_ur and ir_pi in same column
X2sim=squeeze([IR(:,4,INDrnd); IR(:,5,INDrnd)]);

X=squeeze(prctile(Xsim,50,3));
Xlb=squeeze(prctile(Xsim,17,3));
Xub=squeeze(prctile(Xsim,100-17,3));


%median IR estimates
IRmed=squeeze(prctile(IR,50,3));
X1=IRmed(:,1:2); X1=X1(:);%put ir_ur and ir_pi in same column
X2=IRmed(:,4:5); X2=X2(:);

%IR of instruments
X3=IRmed(:,3);
X4=IRmed(:,6);

%To use for color coding of sacatterplot
cc1 = linspace(1,10,H);
cc = linspace(1,10,H);

%Matrices to store OPP distribution
B=nan(TsimEY,TsimX,size(Ey,1));
BIR=nan(TsimX,size(Ey,1));
BEY=nan(TsimEY,size(Ey,1));

% Create a matrix to store counterfactuals
bTablematrix = table2array(bTable);

bTable_output = repmat(bTablematrix, [1, 1, TsimEY*TsimX]);

Ey_gap_adj = zeros(size(Ey_gap, 1), size(Ey_gap, 2), size(Ey_gap, 3));

cd("Counterfactual_2019/");

run("Iteration_SEP_2019");

cd("..")

%% Report results:
bTable_output_SEP = bTable_output;
cd("Counterfactual_2019/");
save('bTable_output_SEP_2019.mat', 'bTable_output_SEP', 'yearn','OPP1','yearn_tp');
Counterfactual_inflation_rate_plot_SEP_2019

toc;