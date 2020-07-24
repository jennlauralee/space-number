function NLL = get_gaussian_oneiter_NLL(X, StartingPoint, mu_resp, conf_resp, params, nSamples)

for i_trial = 1:length(X)
    S = X{i_trial};
    startingpoint = StartingPoint(i_trial);

    [mu_hat_samples, conf_hat_samples] = func_iter_avg_gaussian_oneiter(params, S, startingpoint, nSamples);

%% Fit one Gaussian
    try
        prob(i_trial) = mvnpdf([mu_resp(i_trial),log(conf_resp(i_trial))],[mean(mu_hat_samples) mean(log(conf_hat_samples))], cov(mu_hat_samples, log(conf_hat_samples)));
    catch
        keyboard
    end
end

%% Return the NLL for all trials
prob(prob<exp(-10)) = exp(-10);
NLL = sum(-log(prob));