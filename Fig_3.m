%%
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
% based on non-zero dominate state (1st,3rd,9th and 14th state)
re_occupancy = fractional_occupancy(:,[1,3,9,14]);
re_mean = mean_life(:,[1,3,9,14]);

%% transition calculation

strans_subject_level = compute_subject_level_transition_probabilities(model);

% for rest
transition_subject_wise

for i=1:size(strans_subject_level,2)
   trans_individual(i,:) = reshape(strans_subject_level{1,i}([1,3,9,14],[1,3,9,14])',1,16);
end
figure(1)
figure('Color', 'w');
imagesc(reshape(mean(trans_individual),4,4)')
box off
axis image off



