function [Sc,Lambda,alpha,beta2,gamma2, sig2] = MCMC_AnaFac_Quad_Sparse3_multiregime(Sy,K,a,b,Isnap,Nrun,opt,Ini)
% [Sc,Lambda,alpha2,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3_multiregime(Sy,K,a,b,I,Nrun,opt,Ini)
% Gibbs sampling in the hierachical model (Factor Analysis model) with sparsity enforcement on factors  based on p = 1,...,P measurements
%
% Y(m,i,p) = Lambda(m,k,p)*diag(alpha(k,p))*C(k,i,p) + N(m,i,p), m = 1,...,M, i = 1,...,I
% N(m,i,p) ~ Normal(0,beta2(p)*diag(sig2(m)))
% C(k,i,p) ~ Normal(0,gamma2(p)), k = 1,...,K (K <= M)
% solved from quadratic observation Sy = \sum_i{Y(:,i)*Y(:,i)'}/I
% with
% Lambda(m,k) ~ Normal(0,1/K) (such that E{Lambda*Lambda'} = I)
% sig2(m) ~ InvGamma(a.sig2(m),b.sig2(m)) is the common vector of M noise variances, for the P measurements
% beta2(p) ~ InvGamma(a.beta2(p),b.beta2(p))is the noise relative intensity of the p-th measurement
% alpha(k,p) ~ Exp(a.alpha(k,p))
% gamma2(p) ~ InvGamma(a.gamma2(p),b.gamma2(p))
% which enforces a sparse prior on the raws of C (equivalently on the columns of Lambda).
% opt.noise = 'homo' (homoskedastic): the same sig2 is estimated for all sensors m = 1,..,M 
% opt.noise = 'hetero' (heteroskedastic): a different sig2 is estimated for each sensor m  
%
% Nrun = number of iterations of the Markov chain
% Ini = vector of initial values for Sc, sig2(m), alpha2 and Lambda
%
%If opt.ref = 1, then it is assumed that the 1st measurment is a noise reference : Y(:,:,1)=N(:,:,1)
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: 20/01/2018
%%%%%%%%%%%%%%%%%%%%%%%

M = size(Sy,1);
MK = M*K;
P = size(Sy,3);

if opt.ref == 1 %using a noise reference
	indx=1;
else
	indx=0;
end

if numel(Isnap)==1
	Isnap=Isnap*ones(P,1);
end

for p=1:P
    Sy(:,:,p) = Sy(:,:,p).*Isnap(p);
end

% Check inputs
if ~(strcmp(opt.noise,'homo')||strcmp(opt.noise,'hetero'))
    error('option ''opt.noise'' is wrongly specified!!')
end

if numel(a.alpha) == K
    a.alpha = repmat(a.alpha(:),[1 P-indx]);
    a.alpha = [zeros(K,indx) a.alpha];
elseif numel(a.alpha) == 1
    a.alpha = [zeros(K,indx) a.alpha*ones(K,P-indx)];
elseif numel(a.alpha)==K*(P-indx)
    a.alpha = [zeros(K,indx) a.alpha];
end
    
if numel(a.beta2)==1
   a.beta2=a.beta2*ones(P,1); 
end
if numel(b.beta2)==1
   b.beta2=b.beta2*ones(P,1); 
end

if numel(a.gamma2)==1
   a.gamma2= [zeros(indx);a.gamma2*ones(P-indx,1)];
elseif numel(a.gamma2)==(P-indx)
    a.gamma2=[zeros(indx) a.gamma2];    
else
    disp(['Veuillez donner ' num2str(P-indx) ' valeurs pour a.gamma2']); 
end

if numel(b.gamma2)==1
   b.gamma2= [zeros(indx) ;  b.gamma2*ones(P-indx,1)]; 
elseif numel(b.gamma2)==(P-indx)
    b.gamma2=[zeros(indx) b.gamma2];
else
    disp(['Veuillez donner ' num2str(P-indx) ' valeurs pour b.gamma2']); 
end

if strcmp(opt.noise,'homo')
    a.sig2 = mean(a.sig2); %scalar
    b.sig2 = mean(b.sig2); %scalar
elseif strcmp(opt.noise,'hetero')
    if numel(a.sig2) == 1
        a.sig2 = a.sig2*ones(M,1); %(Mx1)
    end
    if numel(b.sig2) == 1
        b.sig2 = b.sig2*ones(M,1); %(Mx1)
    end
end


% Initialisation
Sc = zeros(Nrun,K,K,P);
Lambda = zeros(Nrun,M,K,P);
alpha = zeros(Nrun,K,P);
beta2 = ones(Nrun,P);
gamma2 = zeros(Nrun,P);
sig2 = zeros(Nrun,M);


if nargin > 7
    Lambda(1,:,:,1+indx:end) = Ini.Lambda; %(NrunxMxKxP)
    if numel(Ini.alpha)~=K*(P-indx); disp(['Veuillez donner ' num2str(K) 'x' num2str(P-indx) ' valeurs pour Ini.alpha']); end;
    alpha(1,:,1+indx:end) = Ini.alpha; %(NrunxKxP)
    sig2(1,:) = Ini.sig2; %(NrunxM)
    if numel(Ini.gamma2)~=P-indx; disp(['Veuillez donner ' num2str(P-indx) ' valeurs pour Ini.gamma2']); end;
    gamma2(1,1+indx:end) = Ini.gamma2; %(NrunxP)
else
	Lambda(1,:,:,1+indx:end) = (randn(M,K,P-indx) + 1i*randn(M,K,P-indx))/sqrt(2); %(NrunxMxKxP)
    alpha(1,:,1+indx:end) = 1./a.alpha(:,1+indx:end);
    sig2(1,:) = b.sig2(:)'./a.sig2(:)'; 
    gamma2(1,1+indx:end) = b.gamma2(1+indx:end)./a.gamma2(1+indx:end);
end

%Initalise with noise reference
%if opt.ref==1
%    sig2(1,:)=real(diag(squeeze(Sy(:,:,1))));
%end

if isfield(opt,'gamma2')
    if opt.gamma2 > 1
        gamma2 = opt.gamma2*ones(Nrun,P);
    end
else
    opt.gamma2 = 0;
end

if strcmp(opt.noise,'hetero')
    a_sig2 = a.sig2 + sum(Isnap);
    a_beta2 = a.beta2 + Isnap(:).*M;
elseif strcmp(opt.noise,'homo')
    a_sig2 = a.sig2 + sum(Isnap)*M;
    a_beta2 = a.beta2 + Isnap(:).*M;
end
a_gamma2 = a.gamma2 + Isnap(:).*K;

Er = zeros(M,1);
for i = 2:Nrun
    for p=1+indx:P
    	% sample in [Sc|rest]     
    	% ===================
    	D=diag((beta2(i-1,p).*sig2(i-1,:)).^-1);
    	LambdaD = squeeze(Lambda(i-1,:,:,p))'*D;
    	Da = diag(alpha(i-1,:,p));
        Cov_inv = Da*LambdaD*squeeze(Lambda(i-1,:,:,p))*Da + eye(K)/gamma2(i-1,p);
        Cov_inv = (Cov_inv + Cov_inv')/2;
        Cov_inv_sqrt = chol(Cov_inv);
        Cov_sqrt = Cov_inv_sqrt\eye(K);
        Cov = Cov_sqrt*Cov_sqrt';
        Wc = Cov_sqrt*Bartlett_Decomp(Isnap(p),eye(K))*Cov_sqrt';
        temp = Cov*Da*LambdaD*Sy(:,:,p)*LambdaD'*Da*Cov' + Wc;
        Sc(i,:,:,p) = (temp + temp')/2;        
    

    
		% sample in [Lambda|rest]
		% =======================
		Omega_inv = kron(Da*conj(squeeze(Sc(i,:,:,p)))*Da,sparse(D)) + K*speye(MK);
		Omega_inv = (Omega_inv + Omega_inv')/2;
		Omega_inv_sqrt = chol(Omega_inv);
		Omega_sqrt = inv(Omega_inv_sqrt);
		Omega = Omega_sqrt*Omega_sqrt';
		Syc = Sy(:,:,p)*D*squeeze(Lambda(i-1,:,:,p))*Da*Cov;
		temp = D*Syc*Da;
		temp = Omega*temp(:) + Omega_sqrt*(randn(MK,1) + 1i*randn(MK,1))/sqrt(2);
		Lambda(i,:,:,p) = reshape(temp,M,K);
		
		% sample in [beta2|rest]
		% ======================
		Z = eye(M) - squeeze(Lambda(i,:,:,p))*Da*Cov*Da*squeeze(Lambda(i,:,:,p))'*D;
		for k = 1:M
		    Er(k,p) = real(Z(k,:)*Sy(:,:,p)*Z(k,:)' + squeeze(Lambda(i,k,:,p)).'*Da*Wc*Da*conj(squeeze(Lambda(i,k,:,p))));
        end
	    b_beta2 = b.beta2(p) + sig2(i-1,:).^-1*Er(:,p);
        %mean(sig2(i-1,:))
	    beta2(i,p) = 1/rand_gamma(1,1,1/b_beta2,a_beta2(p));

        %beta2(i,p)=1;
		
		
		% sample in [gamma2|rest]
    	% =======================
		if opt.gamma2 == 0
		    b_gamma2 = b.gamma2(p);
		    for k = 1:K
		        b_gamma2 = b_gamma2 + real(Sc(i,k,k,p));
		    end
		    gamma2(i,p) = 1/rand_gamma(1,1,1/b_gamma2,a_gamma2(p));
		end
		
		% sample in [alpha|rest]
		% =======================
    	D=diag((beta2(i,p).*sig2(i-1,:)).^-1);
		LambdaD = squeeze(Lambda(i,:,:,p))'*D;
		for k = 1:K
		    var_alpha = .5/real(LambdaD(k,:)*squeeze(Lambda(i,:,k,p)).'*squeeze(Sc(i,k,k,p)));
		    temp  = Syc(:,k) + alpha(i-1,k,p)*squeeze(Sc(i,k,k,p))*squeeze(Lambda(i,:,k,p)).' - squeeze(Lambda(i,:,:,p))*Da*squeeze(Sc(i,:,k,p)).';
		    temp = LambdaD(k,:)*temp;
		    mu_alpha = var_alpha*(2*real(temp) - a.alpha(k,p));
		    %temp = mu_alpha + sqrt(var_alpha)*randn;
		    alpha(i,k,p)=mu_alpha+sqrt(var_alpha)*trandn(-mu_alpha/sqrt(var_alpha),Inf); %Gaussienne non tronquée entre 0 et +Inf
		    %if temp > 0
		    %    alpha(i,k,p) = temp;
		    %else
		    %    alpha(i,k,p) = alpha(i-1,k,p);
		    %end
		end		
    end  
    
    % sample in [sig2|rest]
	% =======================
    
    for p=1:P
        Da=diag(alpha(i,:,p));
        Z = eye(M) - squeeze(Lambda(i,:,:,p))*Da*Cov*Da*squeeze(Lambda(i,:,:,p))'*D;
		for k = 1:M
		    Er(k,p) = real(Z(k,:)*Sy(:,:,p)*Z(k,:)' + squeeze(Lambda(i,k,:,p)).'*Da*Wc*Da*conj(squeeze(Lambda(i,k,:,p))));
        end       
    end
	if strcmp(opt.noise,'hetero')
	    for k = 1:M
	        b_sig2 = b.sig2(k) + Er(k,:)*beta2(i,:)'.^-1;
	        sig2(i,k) = 1/rand_gamma(1,1,1/b_sig2,a_sig2(k));
	    end
	elseif strcmp(opt.noise,'homo')
	    b_sig2 = b.sig2 + sum(sum(beta2(i,:).^-1.*Er(:,:)));
	    sig2(i,:) = ones(1,M)/rand_gamma(1,1,1/b_sig2,a_sig2);
	end 
    
end


for p=1:P
    Sc(:,:,:,p) = Sc(:,:,:,p)./Isnap(p);
end



