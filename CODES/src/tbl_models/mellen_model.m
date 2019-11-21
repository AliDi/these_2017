function [S]=mellen_model(theta,dx,dy)
% [S]=mellen_model([Uc alphax alphay f],dx,dy)
% S = exp(-2*pi*f/uc*(sqrt((alphax*abs(dx)).^2+(alphay*abs(dy)).^2)-1i*dx));
%Calculate the cross-spectral matrix according to the Mellen model
    uc = theta(1);
    alphax = theta(2);
    alphay = theta(3);
    f = theta(4);
	S = exp(-2*pi*f/uc*(sqrt((alphax*abs(dx)).^2+(alphay*abs(dy)).^2)-1i*dx));
    S=(S+S')/2;
end
