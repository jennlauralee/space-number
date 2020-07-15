%% func_iter_avg
% INPUT: X, alpha, sig, maxiter
% OUTPUT: mu_hat, mu_hat_all, niter

% adapted from weiji's "distance_model.m"

% Plot summary statistics after manually inputting parameters

% Get a mu_hat distribution for a single trial
% Fit a mixture of Gaussians to the distribution -- with 3 components --
% get standard mixture Gaussian (EM etc.) don't write own code for that
% This keeps your mu_hat continuous. 
% A few hundred mu_hats; -- pick probability at that distribution of the
% subject's response

% For every trial, subject, parameter, have to create a large fitted
% distribution

function [mu_hat, conf_hat] = func_iter_avg(params, X)

alpha = exp(params(1));
sig = exp(params(2));
convergence_threshold = exp(params(3));


range = 100;
ns = 1e3;

N    = length(X);

mu_hat_0   = rand * range;
ref        = mu_hat_0; % initialize the reference point

iter = 1;
stepsize = inf;
while (stepsize > convergence_threshold)

    % Generative model
    d = ref - X; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
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
    d_s = bsxfun(@times, sign(d), exp(logabsd_s)); % Exponentiate to put in actual d space-- add back in the signs
    refminusmu_s   = mean(d_s,2); % Take mean across lines; vector of ns by 1 -- strive for refminusmu_s to be zero
    %figure
    %hist(refminusmu_s); % Should look like a normal distribution

    mu_s = ref - refminusmu_s; % Absolute coordinates: posterior samples of mu
    %%% Impose post-hoc interval prior: i.e., discard out-of-band mu samples %%%
    intervalprior = 0<=mu_s & 100>=mu_s; % Find in-bound samples
    mu_s = mu_s(intervalprior);
    refminusmu_s = refminusmu_s(intervalprior); % Keep only those in-bound samples

    refminusmu_mean = mean(refminusmu_s); % Get the mean of the samples to get the posterior mean (over the cursor position relative to true mean)
                                          % If the cursor is to the
                                          % RIGHT of hypothesized mean,
                                          % refminusmu_mean is positive.

    % Update
    mu_hat = ref - alpha * refminusmu_mean; 
                       % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                           % rate before you adjust.

    stepsize = abs(ref-mu_hat); % Update absolute stepsize

    ref = mu_hat; % Make the new reference point the posterior mean

    iter = iter + 1;
end

%%%% CONFIDENCE? %%%%

max_halfconf = min([mu_hat, 100-mu_hat]);
halfconfs = [0:1:max_halfconf];
leftconfs = mu_hat - halfconfs;
rightconfs = mu_hat + halfconfs;

AUC = nan(1,length(halfconfs));
for i_c = 1:length(halfconfs)
    AUC(i_c) = sum(mu_s>leftconfs(i_c) & mu_s<rightconfs(i_c))./length(mu_s);
end

rewardfn = exp(-(halfconfs*2)/20);  

expected_utility = AUC.*rewardfn;

[~,i_conf_hat] = max(expected_utility);
conf_hat = halfconfs(i_conf_hat);

% p_halfconf_ = exp(beta * expected_utility); %softmax readout of utility function
% p_halfconf = p_halfconf_./sum(p_halfconf_);

%%%% %%%% %%%% %%%%

          % Build in the reward function to this
          % refminusmu_mean posterior: and get
          % softmax readout of utility function

          % resp_conf_ = abs(pix2x(mouse_loc_pix[0])-resp_loc)*2   %%
          %         resp_conf is the total range of x captured  by the
          %         rectangle (spanning both sides of the response)
          % pts = 15*np.exp(-resp_conf_/20)    
        
end
    