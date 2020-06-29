%% func_iter_avg
% INPUT: X, alpha, sig, maxiter
% OUTPUT: mu_hat, mu_hat_all, niter

% adapted from weiji's "distance_model.m"

function [mu_hat, mu_hat_all, postsample_std, postsample_bincount, postsample_binedges] = func_iter_avg(X, alpha, sig, maxiter)

range = 100;
ns = 1e5;

N    = length(X);

mu_hat_0   = rand * range;
mu_hat     = mu_hat_0;
mu_hat_all = NaN(1,maxiter);

iter = 1;
while (iter <= maxiter)

    % Generative model
    d = mu_hat - X; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
    logx = log(abs(d)) + sig * randn(1,N); % x is vector of measurement of absolute distances to all the lines
                                           % Constant noise on the log of the distance log normal distributions:
                                           % if you add negative noise, you still end up with a positive x because it's
                                           % exponentiated again
                                           
     % Inference through sampling
    logabsd_s = bsxfun(@plus, logx, sig * randn(ns,N)); % For each line (N lines), draw ns samples
                                                        % Add the sample to the log x (fixed observations)

                                                        % adding noise in log space

                                                        % One column of logabsd_s (histogram) looks like lognormal
                                                        % belief distribution.
    d_s = bsxfun(@times, sign(d), exp(logabsd_s)); % Exponentiate to put in actual d space
    refminusmu_s   = mean(d_s,2); % Take mean across lines; vector of ns by 1
    %figure
    %hist(refminusmu_s); % Should look like a normal distribution

    refminusmu_mean = mean(refminusmu_s); % Get the mean of the samples to get the posterior mean (over the cursor position relative to true mean)
                                          % If the cursor is to the
                                          % RIGHT of hypothesized mean,
                                          % refminusmu_mean is positive.
              % Build in the reward function to this
              % refminusmu_mean posterior: and get
              % softmax readout of utility function
    mu_mean = mu_hat - refminusmu_mean; % Mu_mean is not used (just for reference)
                                        % refminusmu_mean is the error
                                        % signal
    % Update
    mu_hat = mu_hat - alpha * refminusmu_mean; % Make the new reference point the posterior mean
                       % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                           % rate before you adjust. 

    mu_hat_all(iter) = mu_hat;

    iter = iter + 1;
end

mu_hat_all = [mu_hat_0, mu_hat_all];

edges = linspace(-50,50,201);
[postsample_bincount,postsample_binedges] = histcounts(refminusmu_s,edges);
postsample_std = std(refminusmu_s);
end
    