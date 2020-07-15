%% test GMModel

close all; clear;
d = 2;
k = 2;
n = 2000;
[X,z] = mixGaussRnd(d,k,n);
plotClass(X,z);
m = floor(n/2);
X1 = X(:,1:m);
X2 = X(:,(m+1):end);
X = X';

%% GMM fitting
GMModel = fitgmdist(X,2);

%%
figure
y = [zeros(1000,1);ones(1000,1)];
h = gscatter(X(:,1),X(:,2),y);
hold on
gmPDF = @(x1,x2)reshape(get_mgvb_density(model,[x1(:) x2(:)]),size(x1));
g = gca;
axis equal
fcontour(gmPDF,[g.XLim g.YLim])
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off
