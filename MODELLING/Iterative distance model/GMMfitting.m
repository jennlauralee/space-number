%% (Per trial) Gaussian mixture model fitting
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

%% Set BADS parameters
% We set the lower/upper bounds for optimization (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [0 -inf log(0)];
UB = [log(1) inf inf];

% We also set the "plausible" lower/upper bounds used by BADS
PLB = [0 log(0.01) log(1)];
PUB = [log(1) log(0.8) log(10)];

% Initial starting parameters
par0 = rand(size(LB)).*(PUB-PLB) + PLB;

nJDsamples = 500;
run = 1;

tic
[pars_run(run,:), NLL_run(run)] = bads(@(par) get_GMM_NLL(par, X, mu_resp, conf_resp, nJDsamples), par0, LB, UB, PLB, PUB);
toc
