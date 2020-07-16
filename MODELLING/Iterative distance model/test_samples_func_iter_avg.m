%% test number of samples within func_iter_avg
clear
clc

X = [12.6004   49.3806];
params =  [0   -2];

nmeasurements = 1000;

[mu,conf] = func_iter_avg(params, X,nmeasurements);


mean(mu)
std(mu)
mean(conf)
std(conf)