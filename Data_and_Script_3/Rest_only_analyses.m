clear all; 
warning off;
rootpath = '';
addpath(genpath(fullfile(rootpath,'Scripts\GCCA_toolbox_sep21\')))
addpath(genpath(fullfile(rootpath,'Scrits\switchingFC\VB-HMM\Scripts\VB-HMM-GD-rev1')))
addpath(genpath('SPM\spm8_scripts'))
addpath(genpath(fullfile(rootpath,'Scripts/CommunityDetection/third_party/2016_01_16_BCT')));
addpath(genpath(fullfile(rootpath,'Scripts/BayesianHiddenFactorAnalysis/privateUse/combinedVAR')))
addpath(genpath(fullfile(rootpath,'Scripts/BayesianHiddenFactorAnalysis/privateUse/VB-HMMFA-AR-v2')))
addpath(genpath(fullfile(rootpath,'Scripts/BayesianHiddenFactorAnalysis/VB-HMMFA-NBD')))
addpath(genpath(fullfile(rootpath,'Scripts/Utilities')));
addpath(genpath(fullfile(rootpath,'Scripts/CommunityDetection/scripts/')));
addpath(genpath('functions'))
addpath(genpath('group_analysis'))

%% ======================================================================================================;
currpath = fullfile(rootpath,'\group_analysis\'); 
savepath = fullfile(rootpath,'\after_init\');
datapath = fullfile(rootpath,'\data\');
stateSet = {'15'}; 
K = 15; 
rep = 'E1E2'; 
trial = '';
filtered = 0;
sessions = {'encoding_1','encoding_2'};
%% load subject behavioral data
pdata_temp = readtable('behavioral_data.xlsx');
pdata = table2array(pdata_temp);

%% load time series data
load('ts.mat');


nSubj = length(data);   display('>>Start.....');
for subj=1:nSubj
    temp = data{subj}; 
    datan{subj} = preprocess(temp);
end

roi_names2 = roi_names;

%% load the trianed model
load('model_3.mat');
%model = group_model; 
estStatesCell = model.temporal_evolution_of_states;
transitions = model.state_transition_probabilities;

for subj = 1:length(estStatesCell)
    estStatesCell_new{subj} = estStatesCell{subj}(1:end); 
end
estStates = cell2mat(estStatesCell);

%% Compute the occupancy rate of each state and the mean life in subject level for a given replication

[fractional_occupancy, mean_life]  = compute_occupancy_and_mean_life_subject_wise(estStatesCell,15);
dominant_states=model.id_of_dominant_states_group_wise;

% plot state temporal evolution
post_states = estStates;
Labels_subj = reshape(post_states,size(datan{1},2),length(datan))';


figure(7);
imagesc(Labels_subj(:,1:end));
pf2 = gca;
set(pf2,'YDir','normal')
ccc=[0.86,0.07,0.23;0.25,0.87,0.82];
colormap(ccc)
ylabel('subject')
 xlabel('time instance')
 title('time evolution')
box off
% mean amplitude

[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);

ori_state_2 = [0.144870522900470,0.0771525559832585,0.0278026167105020,0.0561754980429101,-0.135895610525742,-1.59181600592263e-06,-0.220367598632271,-0.146461083368813,-0.238407582741785,-0.171545907783524,-0.141003606777476,0.0149692730703691,-0.119618020671709,0.0962225321567622,0.850341988419727,0.518206785463390];

[corrE12_ml, pvalueE12_ml] = corr(groupMean{1,9}', ori_state_2','Type','Pearson')

% correlation


ccc=[219,19,59;63,223,208;24,2,112]./255;

x = pdata(1:end,2);
y = mean_life(:,9); xfit = min(y):0.01:max(y)+0.05; 

figure('Color', 'w');
a = plot(mean_life(:,9),x,'o'); set(a,'MarkerEdgeColor',ccc(1,:),'MarkerFaceColor',ccc(1,:));
[p1,s1] = polyfit(mean_life(:,9),x,1);
[yfit,dy] = polyconf(p1,xfit,s1,'predopt','curve','simopt','on');
line(xfit,yfit,'color',ccc(1,:),'LineWidth',2.5)
hold on

p = []; s = [];
xfill = [xfit, fliplr(xfit)];  
yfill = [yfit+dy, fliplr(yfit-dy)];  
fill(xfill, yfill, ccc(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
box off
xlim([0.1 0.55])

[corrE12_ml, pvalueE12_ml] = corr(mean_life(:,[6 9]), pdata(1:end,2),'Type','Spearman')

% transition 

strans_subject_level = compute_subject_level_transition_probabilities(model);

for i=1:size(strans_subject_level,2)
   trans_individual(i,:) = reshape(strans_subject_level{1,i}([6,9],[6,9])',1,4);
end


[corrE12_pr, pvalueE12_pr] = corr(squeeze(strans_subject_level(1,1,:)), pdata(1:end,2), 'Type', 'Spearman') 

for i=1:4
    [corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,i), pdata(1:end,2), 'Type', 'Spearman');
    cor_all(1,i) = corrE12_pr; cor_all(2,i) = pvalueE12_pr;
end




figure('Color', 'w');
x = pdata(1:end,2);
y = trans_individual(:,4)*100;xfit = min(y):0.01:max(y);   % xfit actually y
color2=[0.25,0.87,0.82];
    b = plot(y,x,'o'); set(b,'MarkerEdgeColor',color2,'MarkerFaceColor',color2);
    [p2,s2] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p2,xfit,s2,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color',color2,'LineWidth',2.5)
    hold on
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)]; 
    yfill = [yfit+dy, fliplr(yfit-dy)];  
    fill(xfill, yfill, color2, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
    hold on 
