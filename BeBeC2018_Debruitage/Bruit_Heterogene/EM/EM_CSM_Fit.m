function [L,sig2,beta2,flag,Sx , diag_all, loglik,diffloglik,sig2_all] = EM_CSM_Fit(Sy,I,K,option,Ini)
% [L,sig2,beta2,flag] = EM_CSM_Fit(Sy,I,K,option,Ini)
%
% EM Identification of the stochastic model
%
% Cy(:,:,p) = L(:,:,p)*L(:,:,p)' + beta2(p)*Diag(sig2)
%
% based on p = 1,...,P measurements, where
% Sy(:,:,p) is a MxM cross-spectral matrix estimate of Cy(:,:,p) on I(p) data samples,
% L(:,:,p) is a MxK factor matrix with K < M,
% beta2(p) is the noise relative intensity of the p-th measurement,
% sig2 is the common vector of M noise variances.
% 
% If option.ref = 1, then it is assumed that 
%
% Cy(:,:,1) = Diag(sig2) 
%
% and estimates of L are returned for p = 2,..,P only.
%
% Ini.L and Ini.Syc are initializations of L and beta2(p).*sig2(p) as obtained from SS_SpecSep.
% option.beta2 = 1: different values of beta2 are estimated, otherwise beta2 = 1 is
% assumed for all p's.
% option.max = maximum number of iterations in the EM algorithm.
% option.rerr = target relative error (on sig2) in the EM algorithm.
% flag.count = achieved number of iterations in the EM algorithm.
% flag.norm = evolution of the relative error in the EM algorithm.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: January 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(Sy,1);
%if K >= M,error('K must be less than the number of variates in Sy !');end
P = size(Sy,3);

I = I(:);
if length(I) == 1
    I = I*ones(P,1);
end

% Estimate beta2 (option_beta2 = 1) or force it to unit (option_beta2 = 0)?
if isfield(option,'beta2')
    if option.beta2 == 1
        option_beta2 = 1;
    else
        option_beta2 = 0;
    end
else
    option_beta2 = 0;
end

% Initialization
Syc = Ini.Syc;
L = Ini.L;

sig2 = zeros(M,1);
beta2 = ones(P,1);

% Formating variables
if isfield(option,'ref')
    if option.ref == 1
        if size(Syc,2) ~= (P-1),error(['Number of columns of Ini.Syc must be ',num2str(P-1),' !']);end
        if size(L,3) ~= (P-1),error(['Number of columns of Ini.L must be ',num2str(P-1),' !']);end
        Syc = [real(diag(squeeze(Sy(:,:,1)))) Syc];
        indx = 1;
    else
        if size(Syc,2) ~= P,error(['Number of columns of Ini.Syc must be ',num2str(P),' !']);end
        if size(L,3) ~= P,error(['Number of columns of Ini.L must be ',num2str(P),' !']);end
        indx = 0;
    end
else
    if size(Syc,2) ~= P,error(['Number of columns of Ini.Syc must be ',num2str(P),' !']);end
    if size(L,3) ~= P,error(['Number of columns of Ini.L must be ',num2str(P),' !']);end
    indx = 0;
end

count = 0;
test = 0;
rnorm = zeros(option.max,1);
Iall = sum(I);
while test == 0
    sig2_old = sig2;
    
    % Estimation of sigma2_k
    sig2 = sum(Syc.*repmat((I./beta2)',M,1),2)/Iall;
    
    % Estimation of beta2
    if option_beta2 == 1
        for p = 2:P
            beta2(p) = mean(Syc(:,p)./sig2);
        end
    end
    
    % Estimation of Lambda
    for p = 1:P-indx
        B = squeeze(L(:,:,p))*squeeze(L(:,:,p))' + beta2(p+indx)*diag(sig2);
        B = squeeze(L(:,:,p))'/B;
        Omega = eye(K) - B*squeeze(L(:,:,p));
        
        SB = squeeze(Sy(:,:,p+indx))*B';
        L(:,:,p) = SB/(Omega + B*SB);
        
        % Estimation of residual spectra Syc
        C = eye(M) - squeeze(L(:,:,p))*B;
        for k = 1:M
            Syc(k,p+indx) = real(C(k,:)*squeeze(Sy(:,:,p+indx))*C(k,:)' + squeeze(L(k,:,p))*Omega*squeeze(L(k,:,p))');
        end
        %Alice : Syc=diag(Sy-L*B*Sy);
    end

    % Stop test
    count = count + 1;
        
    
    %Alice : Sx pour chaque iteartion
    B=L*L' + beta2*diag(sig2);
    B=L'/B;
    Omega=eye(K) - B*L;
    SB=Sy*B';    
    diag_all(:,count)=diag(L*(Omega+B*SB)*L');
    
    %calcul de la vraisemblance
    A=L*L'+diag(sig2);    
    loglik(count)=real(-trace(A^(-1)*Sy)-log(det(A)));
    
    if count >=2
    	diffloglik(count)=norm((loglik(count)-loglik(count-1)))/norm(loglik(count-1));
    end
    
    %rnorm(count) = norm(sig2 - sig2_old)/norm(sig2_old);
    rnorm(count) = norm(sig2 - sig2_old)/norm(sig2_old);
    test = (rnorm(count) <= option.rerr) || (count >= option.max);
    if count <=200
    	test=0;
    end
    
    sig2_all(count)=norm(sig2);
    
end
%loglik=loglik./loglik(1);

flag.count = count;
flag.norm = rnorm(1:count);

if nargin > 4
    Sx = zeros(M,M,P);
    for p = 1:P
        B = squeeze(L(:,:,p))*squeeze(L(:,:,p))' + beta2(p+indx)*diag(sig2);
        B = squeeze(L(:,:,p))'/B;
        Omega = eye(K) - B*squeeze(L(:,:,p));
        SB = squeeze(Sy(:,:,p+indx))*B';
        
        Sx(:,:,p) = squeeze(L(:,:,p))*(Omega + B*SB)*squeeze(L(:,:,p))';
        %Sx2=squeeze(L(:,:,p))*B*squeeze(Sy(:,:,p+indx));
    end
end
