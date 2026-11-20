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
load('Data_and_Script_1(main)\model_1.mat');
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

post_states = estStates;
Labels_subj = reshape(post_states,size(datan{1},2),length(datan))';

B = zeros(size(Labels_subj));
B(Labels_subj == 1)  = 1;
B(Labels_subj == 3)  = 2;
B(Labels_subj == 9)  = 3;
B(Labels_subj == 14) = 4;

figure(7);
imagesc(B(:,1:218));
pf2 = gca;
set(pf2,'YDir','normal')
ccc=[0.09,0.01,0.44;0.25,0.87,0.82;1,0.84,0;0.86,0.07,0.23];
colormap(ccc)
c= colorbar('location','eastoutside')
ylabel(c, 'state')
ylabel('subject')
xlabel('time instance')
title('time evolution')
box off

% get mean amplitude and covariance matrix for each state
[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);


for i=1:K
    groupCov_graph{1,i}=groupCov{1,i}-eye(size(groupCov{1,i}));
    groupCov_graph{1,i}(find(groupCov_graph{1,i} < 0))=0; %exclude negative edge
end

figure('Color', 'w');
clim=[0,0.6];
subplot(2,2,1)
groupCov_graph{1,1}(tril(true(size(groupCov_graph{1,1})),-1)) = NaN; 
groupCov_graph{1,1}(tril(true(size(groupCov_graph{1,1})),0)) = NaN; 
imagesc(groupCov_graph{1,1},'AlphaData', ~isnan(groupCov_graph{1,1}))
box off
axis image off
subplot(2,2,2)
groupCov_graph{1,3}(tril(true(size(groupCov_graph{1,3})),-1)) = NaN; 
groupCov_graph{1,3}(tril(true(size(groupCov_graph{1,3})),0)) = NaN; 
imagesc(groupCov_graph{1,3},'AlphaData', ~isnan(groupCov_graph{1,3}))
box off
axis image off
subplot(2,2,3)
groupCov_graph{1,9}(tril(true(size(groupCov_graph{1,9})),-1)) = NaN; 
groupCov_graph{1,9}(tril(true(size(groupCov_graph{1,9})),0)) = NaN; 
imagesc(groupCov_graph{1,9},'AlphaData', ~isnan(groupCov_graph{1,9}))
box off
axis image off
subplot(2,2,4)
groupCov_graph{1,14}(tril(true(size(groupCov_graph{1,14})),-1)) = NaN; 
groupCov_graph{1,14}(tril(true(size(groupCov_graph{1,14})),0)) = NaN; 
imagesc(groupCov_graph{1,14},'AlphaData', ~isnan(groupCov_graph{1,14}))
box off
axis image off

% mean fc comparison
mean(groupCov_graph{1,9}(triu(true(size(groupCov_graph{1,9})), 1)))

[h,p,ci,stats] = ttest((groupCov_graph{1,9}(triu(true(size(groupCov_graph{1,9})), 1))),(groupCov_graph{1,14}(triu(true(size(groupCov_graph{1,14})), 1))))


% calculate network index, load BCT toolbox
[groupCov, groupMean, subjCov, subjMean]= computeDataCovarianceFromDataUsingOnlyStates(datan,estStatesCell, K);
for i=1:K
    groupCov_graph{1,i}=groupCov{1,i}-eye(size(groupCov{1,i}));
    groupCov_graph{1,i}(find(groupCov_graph{1,i} < 0))=0; %exclude negative edge
end

[mean(efficiency_wei(groupCov_graph{1,1})),mean(efficiency_wei(groupCov_graph{1,3})),mean(efficiency_wei(groupCov_graph{1,9})),mean(efficiency_wei(groupCov_graph{1,14}))]

[mean(clustering_coef_bu(groupCov_graph{1,1})),mean(clustering_coef_bu(groupCov_graph{1,3})),mean(clustering_coef_bu(groupCov_graph{1,9})),mean(clustering_coef_bu(groupCov_graph{1,14}))]

% occupancy and mean life

re_occupancy = fractional_occupancy(:,[1,3,9,14]);
re_mean = mean_life(:,[1,3,9,14]);


figure('Color', 'w');
h = boxplot(re_occupancy,'Whisker', Inf,'Widths',0.3);
hold on
for i = 1:28
    if isprop(h(1), 'LineStyle')
        set(h(i), 'LineStyle', '-');  
    end
end
box off

h=findobj(gca,'Tag','Box');
ccc=[0.86,0.07,0.23;1,0.84,0;0.25,0.87,0.82;0.09,0.01,0.44];
for i= 1:length(h)
patch(get(h(i),'XData'),get(h(i),'YData'),ccc(i,:),'FaceAlpha',0.6);
end

%
figure('Color', 'w');
h = boxplot(re_mean,'Whisker', Inf,'Widths',0.3);
hold on
for i = 1:28
    if isprop(h(1), 'LineStyle')
        set(h(i), 'LineStyle', '-');  
    end
end
box off

h=findobj(gca,'Tag','Box');
ccc=[0.86,0.07,0.23;1,0.84,0;0.25,0.87,0.82;0.09,0.01,0.44];
for i= 1:length(h)
patch(get(h(i),'XData'),get(h(i),'YData'),ccc(i,:),'FaceAlpha',0.6);
end

% fo correlation 

[corrE12_pr, pvalueE12_pr] = corr(re_occupancy, pdata(1:end,2), 'Type', 'Spearman') %Separate for Encoding 1 or 2

[corrE12_ml, pvalueE12_ml] = corr(re_mean, pdata(1:end,2),'Type','Spearman')


% transition

strans_subject_level = compute_subject_level_transition_probabilities(model);

for i=1:size(strans_subject_level,2)
   trans_individual(i,:) = reshape(strans_subject_level{1,i}([1,3,9,14],[1,3,9,14])',1,16);
end
figure(1)
figure('Color', 'w');
imagesc(reshape(mean(trans_individual),4,4)')
box off
axis image off

for i=1:16
    [corrE12_pr, pvalueE12_pr] = corr(trans_individual(:,i), pdata(1:end,2), 'Type', 'Spearman');  % pdata(:,1/2) for dprime/overall_accu
    cor_all(1,i) = corrE12_pr; cor_all(2,i) = pvalueE12_pr;
end

trans_all_p = reshape(cor_all(1,:),4,4)';
mafdr(reshape(trans_all_pp(1:4,4),1,4), 'BHFDR', true) 

mafdr(reshape(trans_all_p([3 7 9:12 15]),1,7), 'BHFDR', true)  % fdr for all s3's transition pathway

mafdr(reshape(trans_all_p([4 8 12 13:16]),1,7), 'BHFDR', true)  % fdr for all s4's transition pathway

%% for rest analyses

load('decoding_post_rest.mat','label_data_thr')

[fractional_occupancy_merge, mean_life_merge_r2]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},15);

fo_merge = fractional_occupancy_merge(:,[1,3,9,14]);
fo_merge_r2 = fo_merge./sum(fo_merge,2); % if use threshold

[corrE12_pr, pvalueE12_pr] = corr(mean_life_merge_r2(:,[1,3,9,14]), pdata(1:end,1), 'Type', 'Spearman')

load('decoding_pre_rest.mat','label_data_thr')

[fractional_occupancy_merge, mean_life_merge_r1]  = compute_occupancy_and_mean_life_subject_wise(label_data_thr{1,1},15);

fo_merge = fractional_occupancy_merge(:,[1,3,9,14]);
fo_merge_r1 = fo_merge./sum(fo_merge,2); % if use threshold

[corrE12_pr, pvalueE12_pr] = corr(mean_life_merge_r1(:,[1,3,9,14]), pdata(1:end,2), 'Type', 'Spearman') 

h_temp(:,[1 3 5 7]) = mean_life_merge_r1(:,[1,3,9,14]);
h_temp(:,[2 4 6 8]) = mean_life_merge_r2(:,[1,3,9,14]);

figure('Color', 'w');
h = boxplot(h_temp,'Whisker', Inf,'Widths',0.3);
hold on;
for i = 1:56
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-');  
end
end

[h,p,ci,stats] = ttest(mean_life_merge_r1(:,9),mean_life_merge_r2(:,9))

h_temp(:,[1 3 5 7]) = fo_merge_r1(:,1:4);
h_temp(:,[2 4 6 8]) = fo_merge_r2(:,1:4);

figure('Color', 'w');
h = boxplot(h_temp,'Whisker', Inf,'Widths',0.3);

hold on;
for i = 1:56
if isprop(h(1), 'LineStyle')
set(h(i), 'LineStyle', '-'); 
end
end

[h,p,ci,stats] = ttest(fo_merge_r1(:,3),fo_merge_r2(:,3))

% transition 


% transition analysis for rest
load('decoding_post_rest.mat','label_data_thr')

for i=1:24
trans_rest_r2(:,:,i) = transition_subject_wise(label_data_thr{1,1}{1,i},[1,3,9,14]);
end

imagesc(reshape(mean(trans_rest_r2,3),4,4))

[corrE12_pr, pvalueE12_pr] = corr(squeeze(trans_rest_r2(2,2,:)), pdata(1:end,2), 'Type', 'Spearman') 


load('decoding_pre_rest.mat','label_data_thr')
for i=1:24
trans_rest_r1(:,:,i) = transition_subject_wise(label_data_thr{1,1}{1,i},[1,3,9,14]);
end

[corrE12_pr, pvalueE12_pr] = corr(squeeze(trans_rest_r1(2,2,:)), pdata(1:end,2), 'Type', 'Spearman') 


% permutation test for network analyses
%% Frobenius norm Delta_i = A2_i - A1_i
for i=1:24
  trans_encoding(:,:,i) = strans_subject_level{1,i}([1,3,9,14],[1,3,9,14]);
end

Delta = trans_rest_r2- trans_encoding;   % nNode x nNode x nSub

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

% Individual-specific association of brain state distributions between encoding and post-encoding rest

subplot(2,2,1)
bar(fractional_occupancy(4,[1 3 9]))
subplot(2,2,3)
bar(fo_merge_r2(4,1:3))
subplot(2,2,2)
bar(fractional_occupancy(2,[1 3 9]))
subplot(2,2,4)
bar(fo_merge_r2(2,1:3))


[coef_act, pval]=corr([fractional_occupancy(:,1);fractional_occupancy(:,3);fractional_occupancy(:,9)],[fo_merge(:,1);fo_merge(:,2);fo_merge(:,3)])

for i=1:5000
    rand_order=randperm(24);
    [coef, pval]=corr([fractional_occupancy(rand_order,1);fractional_occupancy(rand_order,3);fractional_occupancy(rand_order,9)],[fo_merge(:,1);fo_merge(:,2);fo_merge(:,3)]);
    coef_null(i)=coef;
end
length(find(coef_null(:)<coef_act))/5000



 clear p1 s1 yfit dy xfill yfill xfit
 y = [fractional_occupancy(:,1),fractional_occupancy(:,3),fractional_occupancy(:,9)];xfit = min(y):0.01:0.5;
 x= [fo_merge_r2(:,1),fo_merge_r2(:,2),fo_merge_r2(:,3)];
color1=[0.25,0.87,0.82];
     a = plot(y,x,'o'); set(a,'MarkerEdgeColor',color1,'MarkerFaceColor',color1);
    [p1,s1] = polyfit(y,x,1); 
    [yfit,dy] = polyconf(p1,xfit,s1,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color','k','LineWidth',2.5)
    hold on
    %plot(xfit,yfit-dy,'k:')
    %plot(xfit,yfit+dy,'k:')
    p = []; s = [];
    xfill = [xfit, fliplr(xfit)]; 
    yfill = [yfit+dy, fliplr(yfit-dy)]; 
    fill(xfill, yfill, color1, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on 
box off

% entropy analyses

for i=1:24
    z = entropy_state(fractional_occupancy(i,:));
    switch_state_encoding(i,1) = z;
   z = entropy_state(fo_merge_r2(i,:));
    switch_state_r2(i,1) = z;
end


figure;
scatter(repmat(1,1,24), switch_state_encoding,50,[1,0.6,0], 'filled'); 
hold on;
scatter(repmat(1.1,1,24), switch_state_r2,50,[0,0.45,0.74], 'filled'); 

for i = 1:length(repmat(1,1,24))
    plot([1, 1.1], [switch_state_encoding(i), switch_state_r2(i)], 'k--'); 
end

[h,p,ci,stats] = ttest(switch_state_r2,switch_state_encoding)

% DDV evidence analyses

for i=1:24
    Labels_subj_r2(i,:) = label_data_thr{1,1}{1,i};
end

post_states = estStates;
Labels_subj_r2 = reshape(post_states,size(datan{1},2),length(datan))';

A = zeros(size(Labels_subj_r2));
A(Labels_subj_r2 == 1)  = 1;
A(Labels_subj_r2 == 3)  = 2;
A(Labels_subj_r2 == 9)  = 3;
A(Labels_subj_r2 == 14) = 4;

figure(7);
imagesc(A(:,1:170));
pf2 = gca;
set(pf2,'YDir','normal')
ccc=[1,1,1;0.09,0.01,0.44;0.25,0.87,0.82;1,0.84,0;0.86,0.07,0.23];
colormap(ccc)
c= colorbar('location','eastoutside')
ylabel(c, 'state')
ylabel('subject')
xlabel('time instance')
title('time evolution')
box off



% evidence 2 compare mean log value for every above threshold point
load('decoding_post_rest.mat','group_log_likelihood')
load('decoding_post_rest.mat','label_data_thr')
for subj = 1:24
    for ns = 1:length(group_log_likelihood(1).data)
        separate_log_likelihood{1,subj}(ns,:) = zscore(group_log_likelihood(1).data{ns}{subj});
    end
end

for i=1:24
    cnt=1;
      temp_pos = find(label_data_thr{1,1}{1,i}(:) ~=0);
    for j=1:length(temp_pos)
       temp=maxk(separate_log_likelihood{1,i}(:,temp_pos(j)),4);
       all_r2{1,i}(cnt) =  (temp(1)-mean(temp))/(std([temp(1) temp(2) temp(3) temp(4)]));
       cnt = cnt+1;
    end
end

mean_r2 = cellfun(@mean, all_r2);
plot(1:24,mean_r2,'-o')
ylim([0.8 1.2])

% permutation
for i=1:24
    for ns = 1:length(group_log_likelihood(1).data)
        separate_log_likelihood{1,i}(ns,:) = zscore(group_log_likelihood(1).data{ns}{i});
    end
    temp_pos = find(label_data_thr{1,1}{1,i}(:) ~=0);
    for q=1:1000
        separate_log_likelihood{1,i}(1,temp_pos) = separate_log_likelihood{1,i}(1,temp_pos(randperm(length(temp_pos))));
        separate_log_likelihood{1,i}(2,temp_pos) = separate_log_likelihood{1,i}(2,temp_pos(randperm(length(temp_pos))));
        separate_log_likelihood{1,i}(3,temp_pos) = separate_log_likelihood{1,i}(3,temp_pos(randperm(length(temp_pos))));
        separate_log_likelihood{1,i}(4,temp_pos) = separate_log_likelihood{1,i}(4,temp_pos(randperm(length(temp_pos))));
        cnt=1;
        for j=1:length(temp_pos)
            
            temp=maxk(separate_log_likelihood{1,i}(:,temp_pos(j)),4);
            
            all_r22{1,i}(cnt) =  (temp(1)-mean(temp))/std([temp(1) temp(2) temp(3) temp(4)]);
            cnt = cnt+1;
        end
        mean_r2_permu(q,i) = mean(all_r22{1,i}(:));
    end
    i
end

length(find(mean(mean_r2_permu,2) > mean(mean_r2)))/1000

% state 2 for post vs pre rest
clear group_log_likelihood label_data_thr all_r1 all_r2
load('decoding_post_rest.mat','group_log_likelihood')
load('decoding_post_rest.mat','label_data_thr')


for subj = 1:24
    for ns = 1:length(group_log_likelihood(1).data)
        separate_log_likelihood{1,subj}(ns,:) = zscore(group_log_likelihood(1).data{ns}{subj});
    end
end

for i=1:24
    cnt=1;
    temp_pos = find(label_data_thr{1,1}{1,i}(:) ==3);
    for j=1:length(temp_pos)
        temp=maxk(separate_log_likelihood{1,i}(:,temp_pos(j)),4);
        all_r2{1,i}(cnt) =  (temp(1)-mean(temp))/(std([temp(1) temp(2) temp(3) temp(4)]));
        cnt = cnt+1;
    end
end

mean_r2 = cellfun(@mean, all_r2);

clear group_log_likelihood label_data_thr
load('decoding_pre_rest.mat','group_log_likelihood')
load('decoding_pre_rest.mat','label_data_thr')

for subj = 1:24
    for ns = 1:length(group_log_likelihood(1).data)
        separate_log_likelihood{1,subj}(ns,:) = zscore(group_log_likelihood(1).data{ns}{subj});
    end
end

for i=1:24
    cnt=1;
    temp_pos = find(label_data_thr{1,1}{1,i}(:) ==3);
    for j=1:length(temp_pos)
        temp=maxk(separate_log_likelihood{1,i}(:,temp_pos(j)),4);
        all_r1{1,i}(cnt) =  (temp(1)-mean(temp))/(std([temp(1) temp(2) temp(3) temp(4)]));
        cnt = cnt+1;
    end
end

mean_r1 = cellfun(@mean, all_r1);


mean(mean_r1)
mean(mean_r2)





%% plot
figure;
hold on;
data1=mean_r1;
data2=mean_r2;
%% boxplot
rng(1);                    
jit = 0.05;
x1 = 1 + (rand(length(data2),1)-0.5)*jit;
x2 = 2 + (rand(length(data1),1)-0.5)*jit;
for i = 1:length(data1)
plot([x1(i) x2(i)],...
[data2(i) data1(i)],...
'--',...
'Color',[0.7 0.7 0.7],...
'LineWidth',0.8);
end
scatter(x1,data2,40,'filled',...
'MarkerFaceColor',[0.85 0.33 0.10],...
'MarkerEdgeColor','k');
scatter(x2,data1,40,'filled',...
'MarkerFaceColor',[0 0.45 0.74],...
'MarkerEdgeColor','k');
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Group 1','Group 2'})
ylabel('Value')
box off
set(gca,'FontSize',12)

% permu 

for i=1:5000
temp_order = randperm(48);
total = [mean_r1,mean_r2];
temp_a = total(temp_order(1:24));
temp_b = total(temp_order(25:48));
all_rand(i) = mean(temp_a)-mean(temp_b);
end
length(find(all_rand >(mean(mean_r2-mean_r1))))/5000


