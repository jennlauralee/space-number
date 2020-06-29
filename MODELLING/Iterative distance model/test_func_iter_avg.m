clear
clc

% load sample X stim from subject 1

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

alpha = 1;
sig = 0.3;
maxiter = 30;

mu_hat = nan(length(X),1);
mu_hat_all = nan(length(X),maxiter+1);

tic
for i_trial = 1:length(X)
    [model.mu_hat(i_trial,1), model.mu_hat_all(i_trial,:), ...
     model.postsample_std(i_trial,1), model.postsample_bincount(i_trial,:)] = func_iter_avg(X{i_trial}, alpha, sig, maxiter);
end
toc

model.alpha = alpha;
model.sig = sig;
model.maxiter = maxiter;

save('test_s1_iter_avg.mat','stim','model');