%%
clear all; 
warning off;
rootpath = 'E:\BSDS_project\BSDS\BSDS\';
addpath(genpath(fullfile(rootpath,'Scripts\GCCA_toolbox_sep21\')))
addpath(genpath('F:\SPM\spm8_scripts'))
addpath(genpath(fullfile(rootpath,'Scripts/CommunityDetection/scripts/')));
addpath(genpath('E:\BSDS_project\BSDS\BSDS\functions'))
addpath(genpath('E:\BSDS_project\BSDS\BSDS\group_analysis'))

load('resting_1_ts.mat')
load('resting_2_ts.mat')


fo_merge = fo_merge./sum(fo_merge,2); % if use threshold
fo_merge = fractional_occupancy_merge(:,[1,3,9,14]);
fo_merge = fo_merge./sum(fo_merge,2); % if use threshold

[fractional_occupancy_merge, mean_life_merge]  = compute_occupancy_and_mean_life_subject_wise(merged_label_data,15);

[fractional_occupancy_merge, mean_life_merge]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},15);


[corrE12_pr, pvalueE12_pr] = corr(fo_merge, pdata(1:end,2), 'Type', 'Spearman') %Separate for encoding 1 or 2
[corrE12_pr, pvalueE12_pr] = corr(mean_life_merge_r2(:,[1,3,9,14]), pdata(1:end,2), 'Type', 'Spearman') %Separate for encoding 1 or 2


[corrE12_pr, pvalueE12_pr] = corr(merged_occupancy, pdata(1:end,4), 'Type', 'Spearman') %Separate for encoding 1 or 2


figure('Color', 'w');
%h = boxplot([re_occupancy,fractional_occupancy_merge(:,[1,3,9,14])],'Whisker', Inf,'Widths',0.3);
h = boxplot([re_occupancy,fo_merge],'Whisker', Inf,'Widths',0.3);
hold on
for i = 1:56
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-');  % 用 set 函数设置属性
end
end
box off
figure('Color', 'w');
h = boxplot(mean_life_merge(:,[1,3,9,14]),'Whisker', Inf,'Widths',0.3);
hold on
for i = 1:28
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-');  % 用 set 函数设置属性
end
end
box off


figure('Color', 'w');
h = boxplot([re_mean,mean_life_merge(:,[1,3,9,14])],'Whisker', Inf,'Widths',0.3);
hold on
for i = 1:56
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-');  % 用 set 函数设置属性
end
end
box off

for i=1:24
    %z = entropy_state(merged_occupancy_r2(i,:));
    z = entropy_state(fractional_occupancy(i,:));
    switch_state_r1(i,1) = z;
   % z = entropy_state(merged_occupancy(i,:));
   z = entropy_state(fo_merge(i,:));
    switch_state_r2(i,1) = z;
end


figure;
scatter(repmat(1,1,24), switch_state_r1,50,[1,0.6,0], 'filled'); % 第一组红色
hold on;
scatter(repmat(1.1,1,24), switch_state_r2,50,[0,0.45,0.74], 'filled'); % 第二组蓝色

for i = 1:length(repmat(1,1,24))
    plot([1, 1.1], [switch_state_r1(i), switch_state_r2(i)], 'k--'); % 黑色虚线连接对应点
end

[h,p,ci,stats] = ttest(switch_state_r2,switch_state_r1)



% r1 r2 compare


fo_merge = fo_merge./sum(fo_merge,2); % if use threshold
fo_merge = fractional_occupancy_merge(:,[1,3,9,14]);
fo_merge = fo_merge./sum(fo_merge,2); % if use threshold

[fractional_occupancy_merge, mean_life_merge]  = compute_occupancy_and_mean_life_subject_wise(merged_label_data,15);

[fractional_occupancy_merge, mean_life_merge]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},15);

[corrE12_pr, pvalueE12_pr] = corr(fo_merge, pdata(1:end,2), 'Type', 'Spearman') %Separate for encoding 1 or 2
[corrE12_pr, pvalueE12_pr] = corr(mean_life_merge_r2(:,[1,3,9,14]), pdata(1:end,2), 'Type', 'Spearman') %Separate for encoding 1 or 2


figure('Color', 'w');
h = boxplot([mean_life_merge_r1(:,[1,3,9,14]),mean_life_merge_r2(:,[1,3,9,14])],'Whisker', Inf,'Widths',0.3);
% h = boxplot([fo_merge_r1,fo_merge_r2],'Whisker', Inf,'Widths',0.3);
hold on;
for i = 1:56
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-');  % 用 set 函数设置属性
end
end

% transition analysis for rest
for i=1:24
trans_rest(:,:,i) = transition_subject_wise(label_data_thr{1,1}{1,i},[1,3,9,14]);
end
figure(1)
imagesc(reshape(mean(trans_rest_r1,3),4,4),[0 0.4])
figure(2)
imagesc(reshape(mean(trans_rest_r2,3),4,4),[0 0.4])


% permutation test for network analyses
%% Delta_i = A2_i - A1_i
Delta = trans_rest_r2- trans_encoding;   

Delta = trans_rest_r2- trans_rest_r1;   
nSub = size(trans_rest_r2,3);
meanDiff = mean(Delta, 3);
D_obs = norm(meanDiff, 'fro');

nperm = 5000;                   
D_perm = zeros(nperm,1);
for p = 1:nperm
    signs = (rand(nSub,1) > 0.5)*2 - 1;               
    meanDiff_perm = mean(Delta .* reshape(signs,1,1,[]), 3);
    D_obs_perm(p) = norm(meanDiff_perm, 'fro');
end

length(D_obs_perm(find(D_obs_perm(:)>D_obs)))/5000


for i=1:size(strans_subject_level,2)
   trans_encoding(:,:,i) = (strans_subject_level{1,i}([1,3,9,14],[1,3,9,14]));
end

