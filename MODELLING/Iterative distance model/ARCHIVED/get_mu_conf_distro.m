%% get_mu_conf_distro

%%

function joint_dist = get_mu_conf_distro(params, X_alltrials, nSamples)

joint_dist = nan(nSamples,2);

for i_trial = 1:length(X_alltrials)
    for i_sample = 1:nSamples
        [joint_dist(i_sample, 1), joint_dist(i_sample,2)] = func_iter_avg(params,X_alltrials);
    end
end

end