function NLL = get_GMM_NLL(params, X, mu_resp, conf_resp, nMeasurements)
%% For each trial, sample the joint distribution from func_iter_avg for the given parameters
prob = nan(length(X),1);
for i_trial = 1:length(X)
    x = X{i_trial};
    [mu_hat_samples, conf_hat_samples] = func_iter_avg(params, x, nMeasurements);
    
%% Fit a Gaussian mixture model (GMM) to the joint distribution using mixGaussVb
    try
        prob(i_trial) = mvnpdf([mu_resp(i_trial),log(conf_resp(i_trial))],[mean(mu_hat_samples) mean(log(conf_hat_samples))], cov(mu_hat_samples, log(conf_hat_samples)));
    catch
        keyboard
    end
    %k(i_trial) = max(z);
    %figure
    %plotClass([mu_hat_samples; conf_hat_samples], z);
        
%% Evaluate the probability of the GMM at the subject's response for all trials
    %prob(i_trial) = get_mgvb_density(model,[mu_resp(i_trial), conf_resp(i_trial)]);
    
% %% visualize 
%     figure
%     scatter(mu_hat_samples, conf_hat_samples);
%     hold on
%     scatter(mu_resp(i_trial), conf_resp(i_trial),100, 'r', 'filled');
%     title(['NLL = ' num2str(-log(prob(i_trial)))])
%     
%     ca = gca;
%     xs = linspace(ca.XLim(1), ca.XLim(2),1000); 
%     ys = linspace(ca.YLim(1), ca.YLim(2),1000);
%     [Xs,Ys] = meshgrid(xs,ys);
%     x_hat = [Xs(:) Ys(:)];
%     [probplot, cov_mat] = get_mgvb_density(model,x_hat);
%     probplot = reshape(probplot,length(ys),length(xs));
%     
%     hold on
%     contour(Xs,Ys,probplot,'k', 'LineWidth', 2)
end

%% Return the NLL for all trials
prob(prob<exp(-10)) = exp(-10);

NLL = sum(-log(prob));

end