function x = Chi2(nu,m,n)
% returns a m-times-n matrix of samples drawn from a Chi2 distribution
% with nu degrees of freedom.
%
%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: 9/11/2014
%%%%%%%%%%%%%%%%%%%%%%

x = squeeze(sum(randn(m,n,nu).^2,3));
