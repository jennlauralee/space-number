close all
clear
clc

addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number'))
%% Which sub, block?
num_block = 1;
sub_id = 1;

%% Which parmas?
load('gaussian_oneiter_k_sub1_Num.mat')
params = model.bestpar;

%params = [log(0.1), log(0.9), log(26)];

%% Condition trials FOR SUB

allcsvdata = readtable('allcsvdata.csv');
i_sub = allcsvdata{:,13}==sub_id & allcsvdata{:,2}==num_block;

mu_resp = allcsvdata{i_sub, 5}; 
conf_resp = allcsvdata{i_sub, 6}; 

allXvec = readtable('allXvec.csv');
i_cut = 2:height(allXvec);
X = cellfun(@str2num, allXvec{i_cut,2}, 'UniformOutput', false); 
X = X(i_sub);

StartingPoint_ = readtable('StartingPoint.csv');
i_cut = 2:height(allXvec);
StartingPoint = StartingPoint_{i_cut,2};
StartingPoint = StartingPoint(i_sub);

stim.maxrange = cellfun(@max, X)-cellfun(@min, X);
stim.std = cellfun(@std, X);
stim.mean = cellfun(@mean, X);

stim.X = X;
stim.StartingPoint = StartingPoint;

stim.num_block = num_block;
stim.sub_id = sub_id;

resp.mu_resp = mu_resp;
resp.conf_resp = conf_resp;

%%
[model.mu_pred, model.conf_pred] = modelpred_gaussian_oneiter(params, stim);
visualize_model_pred(model, resp, stim)

%%
function [mu_pred, conf_pred] = modelpred_gaussian_oneiter(params, stim)

%%

X = stim.X;
nSamples = 500;

mu_hat = nan(1,length(X));
conf_hat = nan(1,length(X));

for i_trial = 1:length(X)
    x = X{i_trial};
    startingpoint = stim.StartingPoint(i_trial);
    [mu_pred(i_trial,:), conf_pred(i_trial,:)] = func_iter_avg_gaussian_oneiter(params,x,startingpoint,nSamples);
end



end

%%
function visualize_model_pred(model, resp, stim)
figd;

subplot(2,2,2)
scatter(stim.mean, mean(model.mu_pred,2))
xlabel('stimulus')
ylabel('mean model estimate')

subplot(2,2,3)
scatter(resp.mu_resp, mean(model.mu_pred,2))
xlabel('response')
ylabel('model estimate')

subplot(2,2,4)
scatter(stim.mean, mean(model.mu_pred,2))
xlabel('stimulus')
ylabel('single-trial model estimate')

figd;
suptitle('Sub 1, Number Block')

subplot(2,2,1)
scatter(stim.mean, resp.mu_resp)
xlabel('stimulus')
ylabel('response')

subplot(2,2,2)
scatter(stim.std, abs(model.mu_pred(:,1) - stim.mean))
xlabel('stim std')
ylabel('absolute error')

subplot(2,2,3)
scatter(stim.mean, abs(model.mu_pred(:,1) - stim.mean))
xlabel('true mean')
ylabel('absolute error')

subplot(2,2,4)
scatter(stim.maxrange, abs(model.mu_pred(:,1) - stim.mean))
xlabel('stim max range')
ylabel('absolute error')


%% plotting summ stats for conf_hat
% figd;
% suptitle('model prediction')
% 
% subplot(2,2,1)
% scatter(stim.mean, model.conf_pred)
% xlabel('true mean')
% ylabel('conf estimate')
% 
% subplot(2,2,2)
% scatter(stim.std, model.conf_pred)
% xlabel('stim std')
% ylabel('conf estimate')
% 
% subplot(2,2,3)
% scatter(stim.maxrange, model.conf_pred)
% xlabel('stim max range')
% ylabel('conf estimate')
end