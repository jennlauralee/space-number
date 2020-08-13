close all
clear
clc

%% which params?
alpha = log(1); %log(0.8);
sig = log(0.2); %log(0.5);
convergence_threshold = log(1);
reward_factor = 1000;

model.pars_run = [alpha, sig, convergence_threshold, reward_factor];

%% Which sub, block?
num_block = 1;
sub_id = 1;

allcsvdata = readtable('allcsvdata.csv');
i_sub = allcsvdata{:,13}==sub_id & allcsvdata{:,2}==num_block;

mu_resp = allcsvdata{i_sub, 5}; 
conf_resp = allcsvdata{i_sub, 6}; 

allXvec = readtable('allXvec.csv');
i_cut = 2:height(allXvec);
X = cellfun(@str2num, allXvec{i_cut,2}, 'UniformOutput', false); 
X = X(i_sub);

stim.maxrange = cellfun(@max, X)-cellfun(@min, X);
stim.std = cellfun(@std, X);
stim.mean = cellfun(@mean, X);

stim.X = X;

stim.num_block = num_block;
stim.sub_id = sub_id;

resp.mu_resp = mu_resp;
resp.conf_resp = conf_resp;


%%
[model.mu_pred, model.conf_pred, model.iter] = modelpred_lognormalsingle(model.pars_run, stim);
visualize_model_pred(model,stim, resp);
% 
% save('modelpred_lognormalsingle.mat', 'model')
% modelpred = [model.mu_pred', model.conf_pred', model.iter'];
% T = array2table(modelpred,'VariableNames', {'mu_pred', 'conf_pred', 'iter_pred'});
% writetable(T,'modelpred_lognormalsingle.csv');
%%
function [mu_pred, conf_pred, iter] = modelpred_lognormalsingle(params, stim)

%%

X = stim.X;
alpha = exp(params(1));
sig = exp(params(2));
convergence_threshold = exp(params(3));
reward_factor = params(4);

mu_hat = nan(1,length(X));
conf_hat = nan(1,length(X));
iter = nan(1,length(X));

for i_trial = 1:length(X)
    x = X{i_trial};
    [mu_pred(i_trial), conf_pred(i_trial), iter(i_trial)] = func_iter_avg_lognormal_single(params,x);
end

end

function visualize_model_pred(model, stim, resp)

figd;

subplot(2,2,1)
scatter(stim.mean, resp.mu_resp)
xlabel('stimulus')
ylabel('response')

subplot(2,2,2)
scatter(stim.mean, model.mu_pred)
xlabel('stimulus')
ylabel('mean model estimate')

subplot(2,2,3)
scatter(resp.mu_resp, model.mu_pred)
xlabel('response')
ylabel('model estimate')

subplot(2,2,4)
scatter(stim.mean, model.mu_pred)
xlabel('stimulus')
ylabel('single-trial model estimate')

    
figd;
suptitle('model prediction')

subplot(2,2,1)
scatter(stim.mean, model.mu_pred)
xlabel('true mean')
ylabel('mean estimate')

subplot(2,2,2)
scatter(stim.std, abs(model.mu_pred - stim.mean'))
xlabel('stim std')
ylabel('absolute error')

subplot(2,2,3)
scatter(stim.mean, abs(model.mu_pred - stim.mean'))
xlabel('true mean')
ylabel('absolute error')

subplot(2,2,4)
scatter(stim.maxrange, abs(model.mu_pred - stim.mean'))
xlabel('stim max range')
ylabel('absolute error')


%% plotting summ stats for conf_hat
figd;
suptitle('model prediction')

subplot(2,2,1)
scatter(stim.mean, model.conf_pred)
xlabel('true mean')
ylabel('conf estimate')

subplot(2,2,2)
scatter(stim.std, model.conf_pred)
xlabel('stim std')
ylabel('conf estimate')

subplot(2,2,3)
scatter(stim.maxrange, model.conf_pred)
xlabel('stim max range')
ylabel('conf estimate')
end