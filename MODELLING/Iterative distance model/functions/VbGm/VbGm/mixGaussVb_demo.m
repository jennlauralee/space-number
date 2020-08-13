%% Variational Bayesian for Gaussian Mixture Model
close all; clear;
d = 2;
k = 2;
n = 10000;
%[X,z] = mixGaussRnd(d,k,n);
true_mu1 = [0 0];
true_mu2 = [3 3];
true_sigma1 = [2 -1; -1 1];
true_sigma2 = [5 1; 1 5];

X1 = mvnrnd(true_mu1, true_sigma1, n/2);
X2 = mvnrnd(true_mu2, true_sigma2, n/2);
X = [X1;X2]';
z = [ones(1,n/2) 2*ones(1,n/2)]';

%scatter(X(:,1),X(:,2))
plotClass(X,z);
m = floor(n/2);
X1 = X(:,1:m);
X2 = X(:,(m+1):end);
axis xy
axis equal
%% VB fitting
[y1, model, L] = mixGaussVb(X,2);
figure;
plotClass(X,y1);
figure;
plot(L)
%% Predict testing data
% [y2, R] = mixGaussVbPred(model,X2);
% figure;
% plotClass(X2,y2);

%% jenn testing
xs = -10:0.2:10; ys = -9:0.2:9;
[Xs,Ys] = meshgrid(xs,ys);
x_hat = [Xs(:) Ys(:)];
[prob, cov_mat] = get_mgvb_density(model,x_hat);
prob = reshape(prob,length(ys),length(xs));

figure
h = scatter(X(1,:),X(2,:),'.');
hold on
contour(Xs,Ys,prob,'k', 'LineWidth', 5)
axis equal
axis xy
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')


hold on % plot true pdf
true_prob1 = mvnpdf(x_hat,true_mu1,true_sigma1); 
true_prob2 = mvnpdf(x_hat,true_mu2,true_sigma2); 

true_prob1 = reshape(true_prob1,length(ys),length(xs));
true_prob2 = reshape(true_prob2,length(ys),length(xs));

contour(Xs,Ys, true_prob1, 'r', 'Linewidth', 3); hold on
contour(Xs,Ys, true_prob2, 'r', 'Linewidth', 3);

