function [corr] = f_corr(sig1,sig2,dim)

sig1 = sig1 - mean(sig1,dim);
sig2 = sig2 - mean(sig2,dim);

std1 = sum(sig1.^2,dim);
std2 = sum(sig2.^2,dim);
std = sqrt(std1.*std2);
cov = sum((sig1.*sig2),dim);

corr = cov./std; 
