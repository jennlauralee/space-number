%% test get_gaussian_oneiter_NLL

clear
clc

 X = {[1 90 50]};
 StartingPoint = 50;
 mu_resp = 47;
 conf_resp = 0;
 nSamples = 100;
 
 params = [log(0.01), log(0.01), log(1/20), log(20)];
 
 tic
 NLL = get_gaussian_oneiter_NLL(X, StartingPoint, mu_resp, conf_resp, params, nSamples)
toc