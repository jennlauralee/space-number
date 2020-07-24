close all
clear
clc

%% Load subject stim + data
addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number/'))

allcsvdata = readtable('allcsvdata.csv');
allXvec = readtable('allXvec.csv');

% allcsvdata table contents:
headers = allcsvdata.Properties.VariableNames;

% allcsvdata{:,2} ==1 for num blocks / == 0 for space blocks
% allcsvdata{:,13} == Sub_ID

%% Condition trials FOR SUB1

i_testsub1 = allcsvdata{:,13}==1 & allcsvdata{:,2}==1; % Get indices for trials for subject 1, num_block ==1
mu_resp = allcsvdata{i_testsub1, 5}; % Get mu_est for subject 1, num_block ==1
conf_resp = allcsvdata{i_testsub1, 6}; % Get conf_est for subject 1
X = cellfun(@str2num, allXvec{i_testsub1,2}, 'UniformOutput', false); % Get X stimuli

teststim = readtable('testsub1stim.csv');
sub1teststim = teststim{:,2};

stim.maxrange = cellfun(@max, X)-cellfun(@min, X);
stim.std = cellfun(@std, X);
stim.mean = cellfun(@mean, X);

stim.X = X;

%%
for i_trial = 1:length(X)
    stim.maxrange(i_trial) = max(X{i_trial})-min(X{i_trial});
    stim.std(i_trial) = std(X{i_trial});
    stim.mean(i_trial) = mean(X{i_trial});
end
stim.X = X;

%%

alpha = log(0.8); %log(0.8);
k_sig_scale = log(0.5);
beta = log(7);
lapse = log(0.1);%log(0.1);
lambda_sig = 1;
params = [alpha, k_sig_scale, beta, lapse, lambda_sig];

nMeasurements = 5e2;

model.mu_hat = nan(1,length(X));
model.conf_hat = nan(1,length(X));
tic

nSamples = 1;

for i_trial = 1:length(X)
    x = X{i_trial};
    [model.mu_hat(i_trial), model.conf_hat(i_trial), nIter(i_trial), mu_hat_all{i_trial}] = func_iter_avg_gaussian_maxiter(params,x,nSamples);%func_iter_avg_single(params,x);%func_iter_avg(params, x,nMeasurements);
    
    %model.mu_hat(i_trial) = mean(mu_hat);
    %model.conf_hat(i_trial) = mean(conf_hat); %Take measurement means
end

toc

model.alpha = alpha;
model.k_sig_scale = k_sig_scale;

%% plotting summ stats for mu_hat

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