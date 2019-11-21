function [a,b] = Convert_Beta(mu,sig)
% [a,b] = Convert_Beta(mu,sig)
% Convert the mean and standard-deviation of a Beta distribution into its shape
% parameters a and b
% mu = (a-1)/(a+b-2)
% sig^2 = a*b/((a+b)^2*(a+b+1))
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Dinsenmeyer
%%%%%%%%%%%%%%%%%%%%%%%%%
var=sig^2;
a = ((1 - mu) / var - 1 / mu) * mu ^ 2;
b = a * (1 / mu - 1);
