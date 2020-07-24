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
LB = [log(0),... %alpha
      log(0),... %k_sig_scale
      log(0),... %beta
      log(0),... %lapse
      log(0)]; %lapse sig
UB = [inf,... %alpha
      inf,... %k_sig_scale
      inf,... %beta
      log(1),... %lapse 
      inf]; %lapse sig
  
PLB = [log(0.01),... %alpha
      log(0.01),... %k_sig_scale
      log(1),... %beta
      log(0.01),... %lapse 
      log(0.01)]; %lapse sig
  
PUB = [log(1),... %alpha
      log(1),... %k_sig_scale
      log(100),... %beta
      log(1),... %lapse
      log(10)]; %lapse sig     

% Initial starting parameters
par0 = rand(size(LB)).*(PUB-PLB) + PLB;

nSamples = 500;
run = 1;

tic
[pars_run(run,:), NLL_run(run)] = bads(@(par) get_gaussian_NLL(X, mu_resp, conf_resp, par, nSamples), par0, LB, UB, PLB, PUB);
toc
