clear all; close all;

g=[ 1 2 3]';
h=g;
 w=[0.1 ; 0.3 ; 0.5];
D=[0 1 2 ; 1 0 6 ; 2 6 0];
pmax=w'*D*w;

for iterH = 1:1000
        hOldValue = h;
        H = h*h';
        
        H(~logical(eye(3))) = 0;
        h = 1/sqrt(1+w'*H*w)*(D*w/pmax + H*w);
        if norm(h-hOldValue) < 1e-6
            break;
        end
end
%G=pmax*h*h';
G=pmax*(h*h'-H);