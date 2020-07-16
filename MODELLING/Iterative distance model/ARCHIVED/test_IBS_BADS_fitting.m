%% test IBS BADS fitting

%% Load subject stim & data
clear
clc

addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number/'))

allcsvdata = readtable('allcsvdata.csv');

i_testsub1 = allcsvdata{:,13}==1 & allcsvdata{:,2}==1; % Get indices for trials for subject 1, num_block ==1
sub1muresp = allcsvdata{i_testsub1, 5}; % Get mu_est for subject 1, num_block ==1
sub1confresp = allcsvdata{i_testsub1, 6}; % Get conf_est for subject 1

teststim = readtable('testsub1stim.csv');
sub1teststim = teststim{:,2};
X = cellfun(@str2num,sub1teststim,'UniformOutput',false);

for i_trial = 1:length(X)
    stim.maxrange(i_trial) = max(X{i_trial})-min(X{i_trial});
    stim.std(i_trial) = std(X{i_trial});
    stim.mean(i_trial) = mean(X{i_trial});
end
stim.X = X;


X = X(1:10);
sub1muresp = sub1muresp(1:10);

%%
% [mu_hat, conf_hat] = func_iter_avg(X, alpha, sig, maxiter);
%% Maximum-likelihood estimation (MLE)

% We fit the data via maximum-likelihood estimation using Bayesian Adaptive 
% Direct Search (BADS), a particularly effective optimization algorithm.
% If you do not have the BADS toolbox installed, you can freely download it 
% from here: https://github.com/lacerbi/bads

%params = [logalpha, logsig, convergence_threshold];


% We set the lower/upper bounds for optimization (in particular, note that 
% we set a nonzero lower bound for the lapse rate)
LB = [0 -inf log(0)];
UB = [log(1) inf inf];

% We also set the "plausible" lower/upper bounds used by BADS
PLB = [0 log(0.01) log(1)];
PUB = [log(1) log(0.8) log(10)];

% We define the negative log-likelihood function via a call to IBSLIKE
% (IBSLIKE provides a vectorized implementation of IBS for MATLAB; check
% out IBS_BASIC for a bare bone implementation)

% We set 10 reps for the IBS estimator (see Section 4.4 in the paper) -
% note that this is also the default for IBSLIKE

options_ibs.Nreps = 1;

%nllfun_ibs = @(params) ibslike(@func_iter_avg, params, [sub1muresp, sub1confresp], X,options_ibs);
nllfun_ibs = @(params) ibslike(@func_iter_avg, params, sub1muresp, X,options_ibs);

% As a starting point for the optimization, we draw a sample inside the
% plausible box (in practice you should use multiple restarts!)
%theta0 = rand(size(LB)).*(PUB-PLB) + PLB;

theta_test = [log(1), log(0.5), log(5)];

tic
%ibs_basic(@func_iter_avg, theta0, [sub1muresp, sub1confresp], X);

%ibs_basic(@func_iter_avg, theta_test, sub1muresp, X);

%ibslike(@func_iter_avg, theta0, [sub1muresp, sub1confresp], X,options_ibs);
[pars_run(run,:), NLL_run(run)] = bads(@(par) func_iter_avg(par, X), theta0, LB, UB, PLB, PUB, [], OPTIONS);
%theta_ibs = bads(nllfun_ibs,theta0,LB,UB,PLB,PUB);
toc

% Compare with MLE obtained using the analytical log-likelihood expression

%fprintf('Maximum-likelihood estimation with BADS using the exact log-likelihood...\n');
%fprintf('(press a key to continue)\n');
%pause;

%theta_exact = bads(@(x) psycho_nll(x,S,R),theta0,LB,UB,PLB,PUB);
