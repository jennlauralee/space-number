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

function [mu_hat, conf_hat, iter] = func_iter_avg(params, X, nMeasurements)
% should give a vector nMeasurements long for both mu_hat and conf_hat

alpha = exp(params(1));
sig = exp(params(2));
%convergence_threshold = exp(params(3));


range = 100;
ns = 1e2; % Samples for a single measurement

N_stim    = length(X);

mu_hat_0   = rand(nMeasurements,1) * range;
ref        = mu_hat_0; % initialize the reference point

iter = 1;
maxiter = 30;
%stepsize = inf;
while iter <maxiter %(stepsize > convergence_threshold)

    % Generative model
    d = ref - X; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
    logx = log(abs(d)) + sig * randn(nMeasurements,N_stim); % x is vector of measurement of absolute distances to all the lines
                                           % Constant noise on the log of the distance log normal distributions:
                                           % if you add negative noise, you still end up with a positive x because it's
                                           % exponentiated again

     % Inference through sampling 
     % LOG NORMAL IS JUST ONE POSSIBILITY FOR
     % WEBER'S LAW. When you add in those signs etc, we have to do
     % inference over sampling-- average of a bunch of log normally
     % distributed numbers cannot be calculated using an equation
     
     % We can simplify this problem: right now we have N_stim objects. 
     % Instead of LOG NORMAL we use a GAMMA distribution (similar shape),
     % you ALSO impose that the std is proportional to the mean
     % (reparameterization), then it's more convenient; all the distances
     % on one side, you can calculate the gamma distribution for their
     % average. Dist of average of the 3 on the right is again a
     % distribution
     
     % The average of gamma distributions is gamma distributed; you do this
     % separately for the left and the right, but you still need samples to
     % combine them across left and right.
    samplenoise =  sig * randn(nMeasurements,N_stim,ns); % draw noise for every fixed observation and sample
    logabsd_s = bsxfun(@plus,logx,samplenoise); %Add noise to the fixed observations
                                                     % For each line (N lines), draw ns samples
                                                        % Add the sample to the log x (fixed observations)

                                                        % adding noise in log space

                                                        % One column of logabsd_s (histogram) looks like lognormal
                                                        % belief distribution.
                                                        
    d_s = bsxfun(@times, sign(d), exp(logabsd_s)); % Exponentiate to put in actual d space-- add back in the signs
    
    refminusmu_s   = squeeze(mean(d_s,2)); % Take mean across lines; vector of ns by 1 -- strive for refminusmu_s to be zero
    %figure
    %hist(refminusmu_s); % Should look like a normal distribution

    mu_s = ref - refminusmu_s; % Absolute coordinates: posterior samples of mu
    %%%%%%%%% Impose post-hoc interval prior: i.e., discard out-of-band mu samples %%%
    %intervalprior = 0<=mu_s & 100>=mu_s; % Find IN-BOUND SAMPLES
    %rejectsamples = ~intervalprior;
    
    %refminusmu_s(rejectsamples) =nan; % Replace out-of-bound samples with nans

    refminusmu_mean = mean(refminusmu_s,2);%nanmean(refminusmu_s,2); % Get the mean of the samples to get the posterior mean (over the cursor position relative to true mean)
                                          % If the cursor is to the
                                          % RIGHT of hypothesized mean,
                                          % refminusmu_mean is positive.

    % Update
    mu_hat = ref - alpha * refminusmu_mean; 
                       % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                           % rate before you adjust.

    %stepsize = abs(ref-mu_hat); % Update absolute stepsize

    ref = mu_hat; % Make the new reference point the posterior mean

    iter = iter + 1;
end

mu_hat = min(mu_hat,100);
mu_hat = max(mu_hat,1);

%%%% CONFIDENCE? %%%%
halfconfs = [0:1:50];
N = histc(abs(mu_s - mu_hat)', halfconfs); %Get the # of samples absolute distance away [0:1:50] from mu_hat
AUC = cumsum(N);

rewardfn = exp(-(halfconfs*2)/20);  

expected_utility = bsxfun(@times,AUC,rewardfn');

[~,i_conf_hat] = max(expected_utility);
conf_hat = halfconfs(i_conf_hat);
conf_hat = min(conf_hat,100);
conf_hat = max(conf_hat,1); 

% THIS TRUNCATED EDGE BREAKS DOWN THE GAUSSIAN APPROXIMATION -- need to fit to a Gaussian with peaked edges  
        
end
    