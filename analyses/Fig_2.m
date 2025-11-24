clear all; 
warning off;
rootpath = 'E:\BSDS_project\BSDS\BSDS\';
addpath(genpath(fullfile(rootpath,'Scripts\GCCA_toolbox_sep21\')))
addpath(genpath('F:\SPM\spm8_scripts'))
addpath(genpath(fullfile(rootpath,'Scripts/CommunityDetection/scripts/')));
addpath(genpath('E:\BSDS_project\BSDS\BSDS\functions'))
addpath(genpath('E:\BSDS_project\BSDS\BSDS\group_analysis'))


%% load the trianed model
load(fullfile(modelpath,trial));
%model = group_model; 
estStatesCell = model.temporal_evolution_of_states;
transitions = model.state_transition_probabilities;

for subj = 1:length(estStatesCell)
    estStatesCell_new{subj} = estStatesCell{subj}(1:end); %leaving out the last 6 empty scans
end
estStates = cell2mat(estStatesCell);

[fractional_occupancy, mean_life]  = compute_occupancy_and_mean_life_subject_wise(estStatesCell,15);
dominant_states=model.id_of_dominant_states_group_wise;


%% get state mean and correlation matrics

groupCov,~,subjCov,~]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell,K);

[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);

[groupMean{1,1};groupMean{1,3};groupMean{1,9};groupMean{1,14}]


figure('Color', 'w');
clim=[-0.2,1];
subplot(2,2,1)
imagesc(groupCov{1,1},clim)
subplot(2,2,2)
imagesc(groupCov{1,3},clim)
subplot(2,2,3)
imagesc(groupCov{1,9},clim)
subplot(2,2,4)
imagesc(groupCov{1,14},clim)

% group level network analyses

%graph analysis load BCT toolbox

groupCov_1 = double(abs(groupCov_1) > 0.01);

% remove diagnoal values
for i=1:K
    groupCov_graph{1,i}=groupCov{1,i}-eye(size(groupCov{1,i}));
    groupCov_graph{1,i}(find(groupCov_graph{1,i} < 0))=0;
end

[mean(efficiency_wei(groupCov_graph{1,1})),mean(efficiency_wei(groupCov_graph{1,3})),mean(efficiency_wei(groupCov_graph{1,9})),mean(efficiency_wei(groupCov_graph{1,14}))]

[mean(clustering_coef_bu(groupCov_graph{1,1})),mean(clustering_coef_bu(groupCov_graph{1,3})),mean(clustering_coef_bu(groupCov_graph{1,9})),mean(clustering_coef_bu(groupCov_graph{1,14}))]
