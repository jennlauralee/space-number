%% (Per trial) Gaussian mixture model fitting
close all
clear
clc

%% Load subject stim + data
addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number/'))

% allcsvdata{:,2} ==1 for num blocks / == 0 for space blocks
% allcsvdata{:,13} == Sub_ID

%% Condition trials 

allcsvdata = readtable('allcsvdata.csv');
mu_resp = allcsvdata{:, 5}; % Get mu_est for subject 1, num_block ==1
conf_resp = allcsvdata{:, 6}; % Get conf_est for subject 1

allXvec = readtable('allXvec.csv');
i_test = 2:height(allXvec);
X = cellfun(@str2num, allXvec{i_test,2}, 'UniformOutput', false); 

StartingPoint_ = readtable('StartingPoint.csv');
i_test = 2:height(allXvec);
StartingPoint = StartingPoint_{i_test,2};

stim.maxrange = cellfun(@max, X)-cellfun(@min, X);
stim.std = cellfun(@std, X);
stim.mean = cellfun(@mean, X);

stim.X = X;
stim.StartingPoint = StartingPoint;
%% % function NLL = func_name(X, mu_resp, conf_resp, params, nJDsamples)

%% Load parameters

% alpha = exp(params(1));
% k_sig_scale = exp(params(2));
% beta = exp(params(3));
% lambda = exp(params(4));
% lapse sig

params = [log(0.1),log(8)];%,log(0.1), log(1)];


%% For each trial, sample the joint distribution from func_iter_avg for the given parameters
tic
testtrials = [1:10]; % Must start at 1
prob = nan(length(testtrials),1);
nSamples = 1e3;
for i_trial = testtrials
    i_trial
    S = X{i_trial};
    startingpoint = StartingPoint(i_trial);
    [mu_hat_samples, conf_hat_samples] = func_iter_avg_gaussian_oneiter(params, S, startingpoint, nSamples);

%% Fit a Gaussian for mu_hat_samples
%% Fit one Gaussian
    try
        prob(i_trial) = mvnpdf([mu_resp(i_trial),log(conf_resp(i_trial))],[mean(mu_hat_samples) mean(log(conf_hat_samples))], cov(mu_hat_samples, log(conf_hat_samples)));
    catch
        keyboard
    end
%% Fit a Gaussian mixture model (GMM) to the joint distribution using mixGaussVb
%     [z, model, ~] = mixGaussVb([mu_hat_samples'],length(S)+1);
%     k(i_trial) = max(z);
%     figure
%     plotClass([mu_hat_samples'; conf_hat_samples], z);
        
%% Evaluate the probability of the GMM at the subject's response for all trials
%    prob(i_trial) = get_mgvb_density(model,[mu_resp(i_trial), conf_resp(i_trial)]);
    
%% visualize 
    figure
    scatter(mu_hat_samples, log(conf_hat_samples));
    hold on
    scatter(mu_resp(i_trial), log(conf_resp(i_trial)),100, 'r', 'filled');
    title(['NLL = ' num2str(-log(prob(i_trial)))])
    
%     ca = gca;
%     xs = linspace(ca.XLim(1), ca.XLim(2),1000); 
%     ys = linspace(ca.YLim(1), ca.YLim(2),1000);
%     [Xs,Ys] = meshgrid(xs,ys);
%     x_hat = [Xs(:) Ys(:)];
%     %[probplot, cov_mat] = get_mgvb_density(model,x_hat);
%     %probplot = reshape(probplot,length(ys),length(xs));
%     
%     hold on
%     %contour(Xs,Ys,probplot,'k', 'LineWidth', 2)
end
toc
%% Return the NLL for all trials
prob(prob<exp(-10)) = exp(-10);
NLL = sum(-log(prob));
