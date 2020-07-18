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

function [mu_hat, conf_hat, iter, mu_hat_all] = func_iter_avg_gaussian_maxiter(params, S, nSamples)
% should give a vector nMeasurements long for both mu_hat and conf_hat

maxiter = 10;

% MAXITER VERISON VECTORIZES THE SAMPLES, provides mu_hat and conf_hat
% samples

alpha = exp(params(1));
k_sig_scale = exp(params(2));
beta = exp(params(3));
lambda = 0;
lapse_sig = 1;
%lambda = exp(params(4));
%lapse_sig = exp(params(5));

range = 100;

N    = length(S);
sig = nan(nSamples,N);
range_s = max(S)-min(S);

mu_hat_0   = rand(nSamples,1) * range;
ref        = mu_hat_0; % initialize the reference point

iter = 1;
delta = inf;

mu_hat_all = nan(nSamples,maxiter);
mu_hat_all(:,1) = mu_hat_0;

S = repmat(S,nSamples,1);
S(rand(nSamples,N)<lambda) = nan;

N    = sum(~isnan(S),2);

while (iter<maxiter)
    %sig(find(N==0),1) = lapse_sig;%ones(sum(N==0),size(sig,2)); % Deal with the case where every stimulus line is dropped...
    
    % Generative model
    d = ref - S; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
    sig(find(N~=0),:) = k_sig_scale*abs(d(find(N~=0),:)); %sigma = scaling constant x absolute distance    
    % Inference
    mu_hat = nanmean(S,2) + randn(nSamples,1).*1/N*sqrt(nansum(sig.^2,2));
    mu_hat(find(N==0)) = mu_hat_0(find(N==0)); %For all samples where everything was dropped, set mu_hat to mu_hat_0...

    % Update
    delta = alpha*(mu_hat-ref); 
    ref = ref + delta;                   % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                           % rate before you adjust.
    mu_hat_all(:,iter) = ref;
    mu_hat = ref; % Store last mu_hat
    iter = iter + 1;
end

%% confidence

post_std = nan(nSamples,1);
post_std = (1./N).*sqrt(nansum(sig.^2,2));
post_std(find(N==0)) = lapse_sig;

halfconfs = [0.1:0.2:50];

combined_vector = combvec(halfconfs,post_std');
AUC = normcdf(combined_vector(1,:),0,combined_vector(2,:)) - normcdf(-combined_vector(1,:),0,combined_vector(2,:));
AUC = reshape(AUC, length(halfconfs), length(post_std));
%AUC = normcdf(halfconfs,0,post_std) - normcdf(-halfconfs,0,post_std);
rewardfn = exp(-(halfconfs*2)/20); 

expected_utility = bsxfun(@times,AUC,rewardfn');

p = exp(beta.*expected_utility); %beta is inverse temperature; higher = lower noise
p = p./sum(p);
for i_sample = 1:nSamples
    try
    conf_hat(i_sample) = randsample(halfconfs,1,'true',p(:,i_sample));
    catch
        keyboard
    end
end

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
    