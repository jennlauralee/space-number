%% (Per trial) Gaussian mixture model fitting using func_iter_avg_SINGLE
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

%% % function NLL = func_name(X, mu_resp, conf_resp, params, nJDsamples)

%% Load parameters
params = [log(1), log(0.5),log(8)];

%% Specifications
nMeasurements = 1000;

%% For each trial, sample the joint distribution from func_iter_avg for the given parameters
tic
testtrials = [1:3]; % Must start at 1
prob = nan(length(testtrials),1);
for i_trial = testtrials
    i_trial
    x = X{i_trial};
    for i_measurement = 1:nMeasurements
        [mu_hat_samples(i_measurement), conf_hat_samples(i_measurement), n_iter(i_measurement)] = func_iter_avg_single(params,x);%func_iter_avg(params, x, nMeasurements);
    end
%% Fit one Gaussian
    try
        prob(i_trial) = mvnpdf([mu_resp(1),log(conf_resp(1))],[mean(mu_hat_samples) mean(log(conf_hat_samples))], cov(mu_hat_samples, log(conf_hat_samples)));
    catch
        keyboard
    end
%% Fit a Gaussian mixture model (GMM) to the joint distribution using mixGaussVb
    %[~, model, ~] = mixGaussVb([mu_hat_samples; conf_hat_samples],1);
    %k(i_trial) = max(z);
    %figure
    %plotClass([mu_hat_samples; conf_hat_samples], z);
        
%% Evaluate the probability of the GMM at the subject's response for all trials
    %prob(i_trial) = get_mgvb_density(model,[mu_resp(i_trial), conf_resp(i_trial)]);
    
%% visualize 
    figure
    scatter(mu_hat_samples, log(conf_hat_samples));
    hold on
    scatter(mu_resp(i_trial), log(conf_resp(i_trial)),100, 'r', 'filled');
    title(['NLL = ' num2str(-log(prob(i_trial)))])
    
    ca = gca;
    xs = linspace(ca.XLim(1), ca.XLim(2),1000); 
    ys = linspace(ca.YLim(1), ca.YLim(2),1000);
    [Xs,Ys] = meshgrid(xs,ys);
    x_hat = [Xs(:) Ys(:)];
    [probplot] = mvnpdf(x_hat,[mean(mu_hat_samples) mean(log(conf_hat_samples))], cov(mu_hat_samples, log(conf_hat_samples)));
    probplot = reshape(probplot,length(ys),length(xs));
    
    hold on
    contour(Xs,Ys,probplot, 'LineWidth', 2)
    colorbar
end
toc
%% Return the NLL for all trials
NLL = sum(-log(prob));
