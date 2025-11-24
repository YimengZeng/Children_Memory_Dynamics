function strans_subject_level = compute_subject_level_transition_probabilities(group_model)
net = group_model.net;

nSubjs = length(net.params.Ybar_old);
nStates = size(net.logOutProbs{1}, 1);
ua = ones(1,nStates)*(1/nStates);
wa0 = zeros(nStates,nStates);
stran = {};
for ns = 1:nSubjs
    phthtpgV1T1 = net.hidden.Qnss{ns}; 
    wa = wa0 + sum(phthtpgV1T1,3); 
    Wa = wa + repmat(ua,[nStates 1]); 
    stran{ns} = exp(  psi(Wa) - repmat( psi(sum(Wa,2)) ,[1 nStates])  );
    sum_stran = sum(stran{ns},2);
    for state=1:nStates
      stran{ns}(state, :) = stran{ns}(state,:)/sum_stran(state);
    end
end

strans_subject_level = stran;
