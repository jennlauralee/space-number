%% get_mu_conf_distro

%% set params
theta_test = [log(1), log(0.5), log(5)];

addpath(genpath('/Users/jennlauralee/GitHub Repos/space-number/'))

teststim = readtable('testsub1stim.csv');
sub1teststim = teststim{:,2};
X = cellfun(@str2num,sub1teststim,'UniformOutput',false);

X = X(1:10);

nSamples = 500;

tic
joint_dist = get_mu_conf_distro(theta_test,X,nSamples);
toc

%%

function joint_dist = get_mu_conf_distro(params, X_alltrials, nSamples)

joint_dist = nan(nSamples,2);

for i_trial = 1:length(X_alltrials)
    for i_sample = 1:nSamples
        [joint_dist(i_sample, 1), joint_dist(i_sample,2)] = func_iter_avg(params,X_alltrials{i_trial});
    end
end
end

%% EM fitting

function EM_fit = fit_EM(

end