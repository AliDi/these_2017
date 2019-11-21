function [S]=corcos_model(theta,dx,dy)
% [S]=corcos_model([Uc alphax alphay f],dx,dy)
% S = exp(-2*pi*f/uc*(alphax*abs(dx)+alphay*abs(dy)-1i*dx));
%Calculate the cross-spectral matrix according to the Corcos model

uc = theta(1);
alphax = theta(2);
alphay = theta(3);
f = theta(4);
S = exp(-2*pi*f/uc*(alphax*abs(dx)+alphay*abs(dy)-1i*dx));
S=(S+S')/2;
end
