%% Func_iter_avg_gaussian_oneiter

% Gaussian noise with a set starting point and a single iteration
% FLEXIBLE noise that scales with distance
% Takes both stimulus and cursor starting position as the reference point
% Alpha (update) is simply 1. Just a one-step update from the starting
% point to the best estimate of the model.

function [mu_hat, conf_hat] = func_iter_avg_gaussian_oneiter(params, S, StartingPoint, nSamples)
% should give a vector nMeasurements long for both mu_hat and conf_hat

k_sig_scale = exp(params(1));
%b_sig = exp(params(2));
beta = exp(params(2));

N    = length(S);

mu_hat_0   = StartingPoint;
ref        = mu_hat_0; % initialize the reference point

S = repmat(S,nSamples,1);

N    = sum(~isnan(S),2);

% Generative model
d = ref - S; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
%sig(find(N~=0),:) = k_sig_scale*abs(d(find(N~=0),:)); %sigma = scaling constant x absolute distance    
sig = k_sig_scale.*abs(d);

% Inference
mu_hat = mean(S,2) + randn(nSamples,1).*1/N*sqrt(mean(sig.^2,2));

% Update
delta = (mu_hat-ref); 
ref = ref + delta;                   % * no alpha (alpha = 1)  % Error signal (refminusmu_mean) should be multiplied by an error
                                       % rate before you adjust.
mu_hat = ref; % Store last mu_hat

%% confidence

post_std = (1./N).*sqrt(sum(sig.^2,2));

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

end
    