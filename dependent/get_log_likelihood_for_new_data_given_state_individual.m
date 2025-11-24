
function logOutProbs = get_log_likelihood_for_new_data_given_state_individual(new_data, bsds_model, pca_model, target_state)

weights = pca_model.weights;
pc_vectors = pca_model.pc_vectors;
scales = sort(weights, 'descend');
rscales = sort(scales(1:size(pc_vectors,2)), 'ascend');
for subj_i=1:length(new_data);
    for pc_i=1:size(pc_vectors,2)
        proj_data(pc_i, :) = rscales(pc_i)*pc_vectors(:,pc_i)'*new_data{subj_i};
    end
    proj_new_data{subj_i} = proj_data;
end

%net = bsds_model.net; % only for group data;take out for individual data
net = bsds_model;
  
Ycell = proj_new_data;
Y = cell2mat(Ycell);
t = target_state;
n = size(Y,2);
p = size(Y,1);

Lm =  net.params.Lm;
psii = net.hparams.psii;
Lcov = net.params.Lcov;
Xm = net.hidden.Xm;
Xcov = net.hidden.Xcov;

% nSubjs = length(Ycell);
for ns = 1:nSubjs
    col = 0; logQns = [];
        col = col + 1;
        kt = size(Lm{t},2);
        LmpsiiLm = Lm{t}'*diag(psii+eps)*Lm{t};
        temp = LmpsiiLm + reshape(reshape(Lcov{t},kt*kt,p)*psii,kt,kt);
        tempXm = Xm{t}(:,1+(ns-1)*n/nSubjs:(ns*n/nSubjs));
        logQns(:,col) = -.5*( +sum(Ycell{ns}.*(diag(psii+eps)*(Ycell{ns}-2*Lm{t}*tempXm)),1)' ...
            +reshape(temp,kt*kt,1)'*reshape(Xcov{t},kt*kt,1) ...
            +sum( tempXm.*(temp*tempXm) ,1)' ...
            +trace(Xcov{t}(2:end,2:end)) ...
            +sum( tempXm(2:end,:).*tempXm(2:end,:) ,1)' ...
            -2*sum(log(diag(chol(Xcov{t}(2:end,2:end)+eps)))) ...
            );
    logOutProbs{ns} = logQns';
end
