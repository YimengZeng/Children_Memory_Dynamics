%% supplementary analyses

% fractional_occupancy and mean life 

domi_state_id = [1,3,9,14];

[fractional_occupancy_merge, mean_life_merge]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},15);
fo_merge = fractional_occupancy_merge(:,domi_state_id);
fo_merge = fo_merge./sum(fo_merge,2); 

[coef_act, pval]=corr([fractional_occupancy(:,domi_state_id(1));fractional_occupancy(:,domi_state_id(2));fractional_occupancy(:,domi_state_id(3))],[fo_merge(:,1);fo_merge(:,2);fo_merge(:,3)])


% temporal sequence correlation between encoding and post-encoding rest
for i=1:5000
    rand_order=randperm(24);
    [coef, pval]=corr([fractional_occupancy(rand_order,1);fractional_occupancy(rand_order,2);fractional_occupancy(rand_order,3);fractional_occupancy(rand_order,4)],[fo_merge(:,1);fo_merge(:,2);fo_merge(:,3);fo_merge(:,4)]);
    coef_null(i)=coef;
end
length(find(coef_null(:)<coef_act))/5000


% occupancy rate entropy
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

% transition in post encoding rest

for i=1:24
trans_rest(:,:,i) = transition_subject_wise(label_data_thr{1,1}{1,i},[1,3,9,14]);
end

Delta = trans_rest_r2- trans_encoding;   % nNode x nNode x nSub

Delta = trans_rest_r2- trans_rest_r1;   % nNode x nNode x nSub
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
