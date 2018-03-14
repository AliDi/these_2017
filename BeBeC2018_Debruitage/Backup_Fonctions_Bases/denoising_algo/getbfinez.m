function b =getbfinez(S)
n=size(S,1);
c=0;
b = zeros(n*(n-1)/2,1);
for i = 1:n-1
    b(c+ [1:(n-i)] )=2*log(abs(S(i+1:end,i)));
    %b(c+ [1:(n-i)] )=S(i+1:end,i);
    c = c+ (n-i);
end

    
