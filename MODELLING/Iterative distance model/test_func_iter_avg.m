clear
clc

% load sample X stim from subject 1

addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number/'))

teststim = readtable('testsub1stim.csv');
teststim = teststim{:,2};
X = cellfun(@str2num,teststim,'UniformOutput',false);


%%
for i_trial = 1:length(X)
    stim.maxrange(i_trial) = max(X{i_trial})-min(X{i_trial});
    stim.std(i_trial) = std(X{i_trial});
    stim.mean(i_trial) = mean(X{i_trial});
end
stim.X = X;

%%

alpha = log(0.8);
sig = log(0.5);
convergence_threshold = log(5);

params = [alpha, sig, convergence_threshold];

mu_hat = nan(length(X),1);

tic

[model.mu_hat, model.conf_hat] = func_iter_avg(params, X);

toc

model.alpha = alpha;
model.sig = sig;
model.convergence_threshold = convergence_threshold;

save('test_s1_iter_avg.mat','stim','model');

%% plotting summ stats for mu_hat

load('test_s1_iter_avg.mat');

figd;
suptitle('model prediction')

subplot(2,2,1)
scatter(stim.mean, model.mu_hat)
xlabel('true mean')
ylabel('mean estimate')

subplot(2,2,2)
scatter(stim.std, abs(model.mu_hat - stim.mean'))
xlabel('stim std')
ylabel('absolute error')

subplot(2,2,3)
scatter(stim.mean, abs(model.mu_hat - stim.mean'))
xlabel('true mean')
ylabel('absolute error')

subplot(2,2,4)
scatter(stim.maxrange, abs(model.mu_hat - stim.mean'))
xlabel('stim max range')
ylabel('absolute error')


%% plotting summ stats for conf_hat
figd;
suptitle('model prediction')

subplot(2,2,1)
scatter(stim.mean, model.conf_hat)
xlabel('true mean')
ylabel('conf estimate')

subplot(2,2,2)
scatter(stim.std, model.conf_hat)
xlabel('stim std')
ylabel('conf estimate')


subplot(2,2,3)
scatter(stim.maxrange, model.conf_hat)
xlabel('stim max range')
ylabel('conf estimate')