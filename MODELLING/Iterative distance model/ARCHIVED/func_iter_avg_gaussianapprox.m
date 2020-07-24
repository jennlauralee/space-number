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

function [mu_hat, conf_hat, iter, mu_hat_all] = func_iter_avg_gaussianapprox(params, S)
% should give a vector nMeasurements long for both mu_hat and conf_hat

alpha = exp(params(1));
k_sig_scale = exp(params(2));
convergence_threshold = exp(params(3));
beta = exp(params(4));
lambda = exp(params(5));

range = 100;

N    = length(S);
range_s = max(S)-min(S);

mu_hat_0   = rand * range;
ref        = mu_hat_0; % initialize the reference point

iter = 1;
delta = inf;

maxiter = 20;
mu_hat_all = nan(1,maxiter);
mu_hat_all(1) = mu_hat_0;

S = S(rand(1,N)>lambda);

N    = length(S);

while (abs(delta) > convergence_threshold) && (iter<maxiter)
    if N == 0
        mu_hat = mu_hat_0;
        sig = randn*100*range_s;
        break
    end
    % Generative model
    d = ref - S; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
    sig = k_sig_scale*abs(d); %sigma = scaling constant x absolute distance    
    % Inference
    mu_hat = mean(S) + randn.*1/N*sqrt(sum(sig.^2));
    
    % Update
    delta = alpha*(mu_hat-ref); 
    ref = ref + delta;                   % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                           % rate before you adjust.
    mu_hat_all(iter) = ref;
    mu_hat = ref; % Store last mu_hat
    iter = iter + 1;
end

%% confidence

post_std = 1/N*sqrt(sum(sig.^2));

halfconfs = [0.1:0.1:50];

AUC = normcdf(halfconfs,0,post_std) - normcdf(-halfconfs,0,post_std);
rewardfn = exp(-(halfconfs*2)/20); 

expected_utility = AUC.*rewardfn;

p = exp(beta*expected_utility); %beta is inverse temperature; higher = lower noise
p = p/sum(p);
conf_hat = randsample(halfconfs,1,'true',p);

%%%% CONFIDENCE? %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max_halfconf = min([mu_hat, 100-mu_hat]);
% halfconfs = [0:1:max_halfconf];
% leftconfs = mu_hat - halfconfs;
% rightconfs = mu_hat + halfconfs;
% 
% AUC = nan(1,length(halfconfs));
% for i_c = 1:length(halfconfs)
%     AUC(i_c) = sum(mu_s>leftconfs(i_c) & mu_s<rightconfs(i_c))./length(mu_s);
% end
% 
% rewardfn = exp(-(halfconfs*2)/20);  
% 
% expected_utility = AUC.*rewardfn;
% 
% [~,i_conf_hat] = max(expected_utility);
% conf_hat = halfconfs(i_conf_hat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    