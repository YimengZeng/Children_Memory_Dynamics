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
pdata_temp = readtable('E:\BSDS_project\code_data_availability\DATA\behavioral_data.xlsx');
pdata = table2array(pdata_temp);

%% load time series data
load('E:\BSDS_project\code_data_availability\time_series\whole_brain\whole_brain_encoding.mat');


nSubj = length(data);   display('>>Start.....');
for subj=1:nSubj
    temp = data{subj}; 
    datan{subj} = preprocess(temp);
end

roi_names2 = roi_names;

%% load the trianed model
load('E:\BSDS_project\code_data_availability\DATA\Data_and_Script_4\model_4.mat');

estStatesCell = model.temporal_evolution_of_states;
transitions = model.state_transition_probabilities;

for subj = 1:length(estStatesCell)
    estStatesCell_new{subj} = estStatesCell{subj}(1:end); 
end
estStates = cell2mat(estStatesCell);

%% Compute the occupancy rate of each state and the mean life in subject level for a given replication

[fractional_occupancy, mean_life]  = compute_occupancy_and_mean_life_subject_wise(estStatesCell,15);
dominant_states=model.id_of_dominant_states_group_wise;

% dominate states: 1,2,3,7,11,15
[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);

[groupMean{1,1};groupMean{1,2};groupMean{1,3};groupMean{1,7};groupMean{1,11};groupMean{1,15}]


% reorder state ID: active-encoding :s1 as s1,s4 as s2; default-mode: s5 as s3; others: s2 as s4, s3 as s5, s6 as s6;

re_occupancy = fractional_occupancy(:,[1,2,3,7,11,15]);
re_mean = mean_life(:,[1,2,3,7,11,15]);

%% Compute correlation between overall memory performance and dominannt states/mean life time
[corrE12_pr, pvalueE12_pr] = corr(re_occupancy(:,4), pdata(1:end,3), 'Type', 'Spearman') 


for i=1:1000
    temp = re_occupancy(randperm(24),4);
    [corrE12_pr, pvalueE12_pr] = corr(temp, pdata(1:end,3), 'Type', 'Spearman');
    all_corr(i) = corrE12_pr;
end

length(find(all_corr>0.44))/1000
%% transition calculation

strans_subject_level = compute_subject_level_transition_probabilities(model);

for i=1:24
    temp = strans_subject_level{1,i}([1,2,3,7,11,15],[1,2,3,7,11,15]);
    temp(eye(size(temp)) == 1) = 0;
    temp = temp ./ sum(temp, 2);
    P_no_self{i} = temp;
    trans_individual(i,:) = reshape(P_no_self{1,i}',1,36);
end
    
    

for i=1:size(strans_subject_level,2)
   trans_individual(i,:) = reshape(strans_subject_level{1,i}([1,2,3,7,11,15],[1,2,3,7,11,15])',1,36);
end



figure('Color', 'w');
imagesc(reshape(mean(trans_individual),3,3)')
box off
axis image off

for i=1:36
    [corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,i), pdata(1:end,3), 'Type', 'Spearman');
    cor_all(1,i) = corrE12_pr; cor_all(2,i) = pvalueE12_pr;
end

% reorder state ID: active-encoding :s1 as s1,s4 as s2; default-mode: s5 as s3; others: s2 as s4, s3 as s5, s6 as s6;


trans_all_p = reshape(cor_all(1,:),6,6)';

[corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,4), pdata(1:end,3), 'Type', 'Spearman')   % s1 to s2
[corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,10), pdata(1:end,3), 'Type', 'Spearman')  % s4 to s2
[corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,22), pdata(1:end,3), 'Type', 'Spearman')  % s2 self transition

for i=1:1000
    temp = trans_individual(randperm(24),10);
    [corrE12_pr, pvalueE12_pr] = corr(temp, pdata(1:end,3), 'Type', 'Spearman');
    all_corr(i) = corrE12_pr;
end

length(find(all_corr>0.4713))/1000


%% trans to s2 （reordered)
figure('Color', 'w');
x = pdata(1:end,3);
y = trans_individual(:,4)*100;xfit = min(y):0.01:max(y); 
color2=[21,96,130]/255;
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
x = pdata(1:end,3);
y = trans_individual(:,10)*100;xfit = min(y):0.01:max(y); 
color3=[160,43,147]/255;
    c = plot(y,x,'o'); set(c,'MarkerEdgeColor',color3,'MarkerFaceColor',color3);
    [p3,s3] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p3,xfit,s3,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color',color3,'LineWidth',2.5)
    hold on
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)];  
    yfill = [yfit+dy, fliplr(yfit-dy)]; 
    fill(xfill, yfill, color3, 'EdgeColor', 'none', 'FaceAlpha', 0.3);  
box off
x = pdata(1:end,3);
y = trans_individual(:,22)*100;xfit = min(y):0.01:max(y); 
color3=[15,158,223]/255;
    c = plot(y,x,'o'); set(c,'MarkerEdgeColor',color3,'MarkerFaceColor',color3);
    [p3,s3] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p3,xfit,s3,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color',color3,'LineWidth',2.5)
    hold on
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)];  
    yfill = [yfit+dy, fliplr(yfit-dy)];  
    fill(xfill, yfill, color3, 'EdgeColor', 'none', 'FaceAlpha', 0.3);  
box off

xlim([0 70])
    
%%  for decoding rest
% reorder state ID: active-encoding :s1 as s1,s4 as s2; default-mode: s5 as s3; others: s2 as s4, s3 as s5, s6 as s6;

load('decoding_post_rest.mat', 'label_data_thr')

[fractional_occupancy_merge, mean_life_merge_r2]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},K);

fo_merge = fractional_occupancy_merge(:,[1,2,3,7,11,15]);
fo_merge_r2 = fo_merge./sum(fo_merge,2); % if use threshold

[corrE12_pr, pvalueE12_pr] = corr(fo_merge_r2(:,5), pdata(1:end,1), 'Type', 'Spearman') 

for i=1:1000
    temp = fo_merge_r2(randperm(24),5);
    [corrE12_pr, pvalueE12_pr] = corr(temp, pdata(1:end,1), 'Type', 'Spearman');
    all_corr(i) = corrE12_pr;
end

length(find(all_corr>0.3191))/1000




for i=1:24
trans_rest_r2(:,:,i) = transition_subject_wise(label_data_thr{1,1}{1,i},[1,2,3,7,11,15]);
end

imagesc(reshape(mean(trans_rest_r2,3),6,6))

[corrE12_pr, pvalueE12_pr] = corr(squeeze(trans_rest_r2(1,5,:)), pdata(1:end,1), 'Type', 'Spearman') 


%% trans post rest
for i=1:1000
    temp = squeeze(trans_rest_r2(1,5,randperm(24)));
    [corrE12_pr, pvalueE12_pr] = corr(temp, pdata(1:end,1), 'Type', 'Spearman');
    all_corr(i) = corrE12_pr;
end

length(find(all_corr>0.44))/1000




%% trans to s2 （reordered)
figure('Color', 'w');
x = pdata(1:end,1);
y = squeeze(trans_rest_r2(1,5,:))*100;xfit = min(y):0.01:max(y); 
color2=[21,96,130]/255;
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
    
    
x = pdata(1:end,3);
y = squeeze(trans_rest_r2(6,5,:))*100;xfit = min(y):0.01:max(y); 
color3=[160,43,147]/255;
    c = plot(y,x,'o'); set(c,'MarkerEdgeColor',color3,'MarkerFaceColor',color3);
    [p3,s3] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p3,xfit,s3,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color',color3,'LineWidth',2.5)
    hold on
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)];  
    yfill = [yfit+dy, fliplr(yfit-dy)];  
    fill(xfill, yfill, color3, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
box off
