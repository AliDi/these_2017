function W = Bartlett_Decomp(n,Sigma)
% W = Bartlett_Decomp(n,Sigma)
% Draw a realization of complex Wishart matrix with mean value n*Sigma and
% degree of freedom n.
% If Sigma is a vector, it implies a diagonal mean value matrix.
% If a scalar m is specified instead, a mean value n*I*m is assumed.
%
% J. Antoni: 15/04/2012

[m1,m2] = size(Sigma);
m = max(m1,m2);
if m == 1
    m = Sigma;
end

% Synthesis of W(m,I)
T = zeros(m,n);
for i = 1:min(n,m)
    T(i,i) = sqrt(sum(randn(2*(n-i+1),1).^2));
end
for i = 2:m
    p = min(i-1,n);
    T(i,1:p) = randn(p,1) + 1i*randn(p,1);
end
W = T*T'/2;

% Synthesis of W(m,Sigma)
if m == min(m1,m2)
    Sigma = (Sigma + Sigma')/2;
    Sigma = chol(Sigma)';
else
    Sigma = diag(sqrt(Sigma));
end
W = Sigma*W*Sigma';
W = (W+W')/2;

