function [W,D] = cwishrnd(sigma,df,D)
% [W,D] = cwishrnd(sigma,df,D)
% Generate Complex wishart matrix with df degree of freedom and covariance matrix sigma.
%
% W = cwishrnd(sigma,df) Generate Complex wishart matrix with df degree of freedom. If sigma is a square matrix, it is used as covariance matrix. If sigma = [a N] with a and N scalars , the covariance matrix is assumed to be a*eye(N).
%
% W = cwishrnd(sigma,df,D) expects D to be the Cholesky factor of sigma. 
%
% [W,D] = cwishrnd(sigma,df) returns D.

narginchk(2,3)
flag=0;
[n,m]=size(sigma);

if sum([n,m])==3
	n=sigma(2);
	a=sigma(1);
	flag=1;
elseif (n~=m || sum([n m])==2)
	error('sigma must be square or a vector of 2 scalars')
end

T=diag(sqrt(chi2rnd(2*(df-(0:n-1)))));
T(itriu(n)) = randn(n*(n-1)/2,1) + 1i*randn(n*(n-1)/2,1);	
W = T'*T/2;

if flag
	W=a*W;
	D=[];
else
	if nargin<3
		D=chol(sigma);
    end
    W=D'*W*D;
end
end


% --------- get indices of upper triangle of p-by-p matrix
function d=itriu(p)

d=ones(p*(p-1)/2,1);
d(1+cumsum(0:p-2))=p+1:-1:3;
d = cumsum(d);
end
