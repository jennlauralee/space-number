% 201910123 Weiji's iterative averaging model, in which each iteration
% consists of measuring distances from a reference point (the current estimate), calculating the
% posterior over the mean relative to the reference point, and choosing the
% posterior mean as the new reference point.

% Could fit # of iterations (or fixed # of iterations and observe errors)
% Alpha, sigma, and maxiter could all become free parameters
% Analyze behaviour of model where if variance is larger, get bigger errors

% Can get confidence from this model:
% How might max iter and sig trade off against each other?

clear; close all

alpha = 1; % Learning rate (responsivity to error)
maxiter = 30;
range = 100;
sig = 0.3;
ns = 1e5;
ntrials = 5;

figure;
set(gcf,'Position',[0 0 560 420])
ax = gca;
ax.NextPlot = 'replaceChildren';

v = VideoWriter('trials.avi');
open(v);

for t = 1:ntrials
    iter = 1;
    
    N = randi(5)+1; % number of points to average
    s = rand(1,N) * range;
    
    mu_hat_0   = rand * range;
    mu_hat     = mu_hat_0;
    mu_hat_all = NaN(1,maxiter);
    
    while (iter <= maxiter)
        
        % Generative model
        d = mu_hat - s; % d is how far off the cursor is (e.g. positive means cursor is RIGHT of the line)
        logx = log(abs(d)) + sig * randn(1,N); % x is vector of measurement of absolute distances to all the lines
                                               % Constant noise on the log of the distance
                                               % log normal distributions:
                                               % if you add negative noise,
                                               % you still end up with a
                                               % positive x because it's
                                               % exponentiated again
        % Inference through sampling
        logabsd_s = bsxfun(@plus, logx, sig * randn(ns,N)); % For each line (N lines), draw ns samples
                                                            % Add the
                                                            % sample to the
                                                            % log x (fixed
                                                            % observations)
                                                            
                                                            % adding noise
                                                            % in log space
        
                                                            % One column of
                                                            % logabsd_s
                                                            % (histogram)
                                                            % looks like
                                                            % lognormal
                                                            % belief
                                                            % distribution.
        d_s = bsxfun(@times, sign(d), exp(logabsd_s)); % Exponentiate to put in actual d space
                                                        % Put the true signs on every sample
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
                  
                  %calc expected utility of different regions; for each
                  %candidate confidence region size, calc expected reward
                  %(points x probability you'll get reward)
                  % probability you get reward = the proportion of samples
                  % falling within confidence region size
                  
        mu_mean = mu_hat - refminusmu_mean; % Mu_mean is not used (just for reference)
                                            % refminusmu_mean is the error
                                            % signal
        
%         if t == 1 && iter == 1
%             figure
%             hist(refminusmu_s)
%         end
        
        % Update
        mu_hat = mu_hat - alpha * refminusmu_mean; % Make the new reference point the posterior mean
                           % * alpha   % Error signal (refminusmu_mean) should be multiplied by an error
                                               % rate before you adjust. 
        
        mu_hat_all(iter) = mu_hat;
        
        iter = iter + 1;
    end
    mu_hat_all = [mu_hat_0, mu_hat_all];
    
    for j=1:maxiter+20
        clf
        axis off
        axis([0 100 -10 10])
        line([0 100], [0 0],'Color','blue','linewidth', 1)
        
        for i=1:N
            line(s(i)*[1 1], [-1 1],'Color','black','linewidth', 1.5)
        end
        
        if j<maxiter
            line(mu_hat_all(j) * [1 1], [-1 1],'Color','red','linewidth', 1.5)
        else
            line(mu_hat_all(maxiter) * [1 1], [-1 1],'Color','red','linewidth', 1.5)
        end
        if j> maxiter + 10
            line(mu_hat_all(maxiter) * [1 1], [-1 1],'Color','red','linewidth', 1.5)
            line(mean(s) * [1 1], [-1 1],'Color','green','linewidth', 1.5)
        end
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    model.mu_hat(t) = mu_hat;
    stim.mean(t) = mean(s);
end
close(v);