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

num_block = 1;
sub_id = 1;

i_testsub = allcsvdata{:,13}==sub_id & allcsvdata{:,2}==num_block; % Get indices for trials for subject 1, num_block ==1
mu_resp = allcsvdata{i_testsub, 5}; % Get mu_est for subject 1, num_block ==1
conf_resp = allcsvdata{i_testsub, 6}; % Get conf_est for subject 1
X = cellfun(@str2num, allXvec{i_testsub,2}, 'UniformOutput', false); % Get X stimuli

stim.maxrange = cellfun(@max, X)-cellfun(@min, X);
stim.std = cellfun(@std, X);
stim.mean = cellfun(@mean, X);

stim.X = X;

%% Set BADS parameters
% We set the lower/upper bounds for optimization (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [log(0),... %alpha
      log(0),... %k_sig_scale
      log(0)];%,... %beta
      %log(0),... %lapse
      %log(0)]; %lapse sig
UB = [inf,... %alpha
      inf,... %k_sig_scale
      inf];%,... %beta
      %log(1),... %lapse 
      %inf]; %lapse sig
  
PLB = [log(0.01),... %alpha
      log(0.01),... %k_sig_scale
      log(1)];%,... %beta
      %log(0.01),... %lapse 
      %log(0.01)]; %lapse sig
  
PUB = [log(1),... %alpha
      log(1),... %k_sig_scale
      log(100)];%,... %beta
      %log(1),... %lapse
      %log(10)]; %lapse sig     

% Initial starting parameters
par0 = rand(size(LB)).*(PUB-PLB) + PLB;

nSamples = 500;

tic
[pars_run, NLL_run] = bads(@(par) get_gaussian_maxiter_NLL(X, mu_resp, conf_resp, par, nSamples), par0, LB, UB, PLB, PUB);
toc

[~, bestrunidx] = max(NLL_run);


%% save output
model.pars_run = pars_run;
model.NLL_run = NLL_run;
model.sub_id = sub_id;
model.num_block = num_block;

if num_block
    block_text = 'Num';
else
    block_text = 'Space';
end

save(['gaussian_maxiter_sub' num2str(sub_id) '_' block_text '.mat'], 'stim','model');