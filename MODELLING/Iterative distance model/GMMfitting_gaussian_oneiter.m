%% (Per trial) Gaussian mixture model fitting
close all
clear
clc

%%
k_sig = log(0.2); %log(0.8);
beta = log(8);

params = [k_sig,beta];

%% Condition trials FOR SUB1

num_block = 1;
sub_id = 1;

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


%% Set BADS parameters
% We set the lower/upper bounds for optimization (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [log(0),... %k_sig_scale
      log(0)];%,... %beta
UB = [inf,... %k_sig_scale
      inf];%,... %beta
  
PLB = [log(0.01),... %k_sig_scale
      log(1)];%,... %beta
  
PUB = [log(1),... %k_sig_scale
      log(100)];%,... %beta    

% Initial starting parameters
par0 = rand(size(LB)).*(PUB-PLB) + PLB;

nSamples = 500;

tic
[pars_run, NLL_run] = bads(@(par) get_gaussian_oneiter_NLL(stim.X, stim.StartingPoint, resp.mu_resp, resp.conf_resp, params, nSamples), par0, LB, UB, PLB, PUB);
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

save(['gaussian_oneiter_sub' num2str(sub_id) '_' block_text '.mat'], 'stim','resp','model');