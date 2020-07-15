%% calculate probability density from mgvb (mixGaussVb) model 

function [prob, cov_mat] = get_mgvb_density(model,x_hat)

D = size(model.m,1);
K = size(model.m,2);
alpha = model.alpha;
alpha_hat = sum(alpha);
m = model.m;
v = model.v;
beta = model.kappa;
L_prefactor = ((v+1-D).*beta)./(1+beta);
M = model.M;

L_k = nan(D,D,K);
for i_k = 1:K
    W_k = inv(M(:,:,i_k));
    L_k(:,:,i_k) = L_prefactor(i_k).*W_k;
end

cov_mat = nan(D,D,K);
prob_k = nan(K,1);
for i_k = 1:K 
    C = inv(L_k(:,:,i_k)); %Covariance matrix
    cov_mat(:,:,i_k) = C; %Store covariance matrix for output
    s = sqrt(diag(C)); %s is the standard deviation vector
    correlation_mat = diag(1./s)*C*diag(1./s);
    x_hat_standardized = bsxfun(@rdivide, x_hat-m(:,i_k)',s'); % bsxfun with @rdivide will do repmats then do ./s
    prob_k(i_k) = alpha(i_k)*mvtpdf(x_hat_standardized, correlation_mat, v(i_k)+1-D); % eq 10.81 
    % To get the value of a normal distribution, you can do (x - mu) / std,
    % then plug in a standard normal (std = 1, mu = 0)
    % This is what's happening in mvtpdf, where it assumes std = 1, mu =0
    % -- so x_hat has to be standardized (line 25)
end

prob = (1/alpha_hat)*sum(prob_k);

end