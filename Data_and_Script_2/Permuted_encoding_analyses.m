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
stateSet = {'15'}; % number of states
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
load('model_2.mat');
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
imagesc(Labels_subj(:,1:218));
pf2 = gca;
set(pf2,'YDir','normal')
%ccc=[219,19,59;255,215,0;63,223,208;24,2,112]./255;
ccc=[219,19,59;63,223,208;24,2,112]./255;
colormap(ccc)
c= colorbar('location','eastoutside')
ylabel(c, 'state')
ylabel('subject')
 xlabel('time instance')
 title('time evolution')
box off
% correlation 
[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);


[groupMean{1,3};groupMean{1,9};groupMean{1,14}]


% correlation

re_occupancy = fractional_occupancy(:,[3,9,14]);

[corrE12_ml, pvalueE12_ml] = corr(re_occupancy(:,[1:3]), pdata(1:end,2),'Type','Spearman')

mafdr([pvalueE12_ml], 'BHFDR', true)  


ccc=[219,19,59;63,223,208;24,2,112]./255;
x = pdata(1:end,2);
y = re_occupancy(:,1)*100; xfit = min(y):0.01:max(y)+0.05; 

figure('Color', 'w');
a = plot(re_occupancy(:,1)*100,x,'o'); set(a,'MarkerEdgeColor',ccc(1,:),'MarkerFaceColor',ccc(1,:));
[p1,s1] = polyfit(re_occupancy(:,1)*100,x,1);
[yfit,dy] = polyconf(p1,xfit,s1,'predopt','curve','simopt','on');
line(xfit,yfit,'color',ccc(1,:),'LineWidth',2.5)
hold on

p = []; s = [];
xfill = [xfit, fliplr(xfit)];  
yfill = [yfit+dy, fliplr(yfit-dy)]; 
fill(xfill, yfill, ccc(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
hold on
xfit = min(re_occupancy(:,2)*100):0.01:max(re_occupancy(:,2)*100)+0.05; 
b = plot(re_occupancy(:,2)*100,x,'o'); set(b,'MarkerEdgeColor',ccc(2,:),'MarkerFaceColor',ccc(2,:));
[p2,s2] = polyfit(re_occupancy(:,2)*100,x,1);
[yfit,dy] = polyconf(p2,xfit,s2,'predopt','curve','simopt','on');
line(xfit,yfit,'color',ccc(2,:),'LineWidth',2.5)
hold on

xfill = [xfit, fliplr(xfit)];  
yfill = [yfit+dy, fliplr(yfit-dy)];  
fill(xfill, yfill, ccc(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);  
xfit = min(re_occupancy(:,3)*100):0.01:max(re_occupancy(:,3)*100)+0.05; 
c = plot(re_occupancy(:,3)*100,x,'o'); set(c,'MarkerEdgeColor',ccc(3,:),'MarkerFaceColor',ccc(3,:));
[p3,s3] = polyfit(re_occupancy(:,3)*100,x,1);
[yfit,dy] = polyconf(p3,xfit,s3,'predopt','curve','simopt','on');
line(xfit,yfit,'color',ccc(3,:),'LineWidth',2.5)
hold on

xfill = [xfit, fliplr(xfit)];  
yfill = [yfit+dy, fliplr(yfit-dy)]; 
fill(xfill, yfill, ccc(3,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3); 


box off
xlim([0.1 0.55])

% transition 

strans_subject_level = compute_subject_level_transition_probabilities(model);

for i=1:size(strans_subject_level,2)
   trans_individual(i,:) = reshape(strans_subject_level{1,i}([3,9,14],[3,9,14])',1,9);
end

reshape(mean(trans_individual),3,3)'

for i=1:9
    [corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,i), pdata(1:end,2), 'Type', 'Spearman');
    cor_all(1,i) = corrE12_pr; cor_all(2,i) = pvalueE12_pr;
end

mafdr(cor_all(2,[1:3 4 7]), 'BHFDR', true)  %correct for all s1's transition
mafdr(cor_all(2,[2 4:6 8]), 'BHFDR', true)  %correct for all s2's transitio
mafdr(cor_all(2,[3 6 7:9]), 'BHFDR', true)   %correct for all s3's transitio

trans_all_p = reshape(cor_all(1,:),3,3)';

% transition plot
%% trans to s3

[corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,8), pdata(1:end,2), 'Type', 'Spearman')
ccc=[219,19,59;63,223,208;24,2,112]./255;

figure('Color', 'w');
x = pdata(1:end,2);
y = trans_individual(:,1)*100;xfit = min(y):0.01:max(y); 
color2=ccc(1,:);
    b = plot(y,x,'o'); set(b,'MarkerEdgeColor',color2,'MarkerFaceColor',color2);
    [p2,s2] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p2,xfit,s2,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color',color2,'LineWidth',2.5)
    hold on
    %plot(xfit,yfit-dy,'k:')
    %plot(xfit,yfit+dy,'k:')
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)]; 
    yfill = [yfit+dy, fliplr(yfit-dy)];  
    fill(xfill, yfill, color2, 'EdgeColor', 'none', 'FaceAlpha', 0.3);  
    hold on 
x = pdata(1:end,2);
y = trans_individual(:,4)*100;xfit = min(y):0.01:max(y); 
color3=ccc(2,:);
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
x = pdata(1:end,2);
y = trans_individual(:,7)*100;xfit = min(y):0.01:max(y); 
color3=ccc(3,:);
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

    xfill = [xfit, fliplr(xfit)];  
    yfill = [yfit+dy, fliplr(yfit-dy)];  
    fill(xfill, yfill, color3, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
box off

xlim([12 56])


