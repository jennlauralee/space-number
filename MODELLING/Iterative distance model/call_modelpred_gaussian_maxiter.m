
[model.mu_pred, model.conf_pred] = modelpred_gaussian_maxiter(model.pars_run, stim);
visualize_model_pred(model,stim);

%%
function [mu_pred, conf_pred] = modelpred_gaussian_maxiter(params, stim)

%%

X = stim.X;
alpha = exp(params(1));
k_sig_scale = exp(params(2));
beta = exp(params(3));

nMeasurements = 5e2;

mu_hat = nan(1,length(X));
conf_hat = nan(1,length(X));
tic

nSamples = 1;

for i_trial = 1:length(X)
    x = X{i_trial};
    [mu_pred(i_trial), conf_pred(i_trial), nIter(i_trial), mu_hat_all{i_trial}] = func_iter_avg_gaussian_maxiter(params,x,nSamples);
end

end

function visualize_model_pred(model, stim)
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