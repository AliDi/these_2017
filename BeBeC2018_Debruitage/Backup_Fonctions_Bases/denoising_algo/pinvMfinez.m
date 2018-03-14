function out=pinvMfinez(n)
M=[];
for i = 1:n
    M= [M;[ zeros(n-i,i-1) ones(n-i,1) eye(n-i,n-i)]];
end
out = pinv(double(M));
%out=M;

