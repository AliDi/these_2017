function [Sc,Lambda,q , l ,beta2,gamma2, sig2,s] = MCMC_AnaFac_Quad_BernGauss_multiregime(Sy,K,a,b,Isnap,Nrun,opt,Ini)
% [Sc,Lambda,q,l,beta2,gamma2 , sig2] = MCMC_AnaFac_Quad_BernGauss_multiregime(Sy,K,a,b,I,Nrun,opt,Ini)
% Gibbs sampling in the hierachical model (Factor Analysis model) with sparsity enforcement on factors  based on p = 1,...,P measurements
%
% Y(m,i,p) = Lambda(m,k,p)*C(k,i,p) + N(m,i,p), m = 1,...,M, i = 1,...,I
% N(m,i,p) ~ Normal(0,beta2(p)*diag(sig2(m)))
% C(k,i,p) ~ BernGauss(l(i,p),gamma2(k,p)), k = 1,...,K (K <= M)
% solved from quadratic observation Sy = \sum_i{Y(:,i)*Y(:,i)'}/I
% with
% Lambda(m,k) ~ Normal(0,1/K) (such that E{Lambda*Lambda'} = I)
%
%Hyperpriors:
% sig2(m) ~ InvGamma(a.sig2(m),b.sig2(m)) is the common vector of M noise variances, for the P measurements
% beta2(p) ~ InvGamma(a.beta2(p),b.beta2(p))is the noise relative intensity of the p-th measurement
% l(i,p) ~ Beta(a.l(p),b.l(p))) sparsity parameter, associated with the Bernoulli process
% gamma2(k,p) ~ InvGamma(a.gamma2(k,p),b.gamma2(k,p))
% which enforces a sparse prior on the raws of C (equivalently on the columns of Lambda).
% opt.noise = 'homo' (homoskedastic): the same sig2 is estimated for all sensors m = 1,..,M
% opt.noise = 'hetero' (heteroskedastic): a different sig2 is estimated for each sensor m
%
% Nrun = number of iterations of the Markov chain
% Ini = vector of initial values for Sc (optional), sig2(m), gamma2 and Lambda
%
%If opt.ref = 1, then it is assumed that the 1st measurment is a noise reference : Y(:,:,1)=N(:,:,1)
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: 20/01/2018,
% Update A. Dinsenmeyer,
% from C. Faure thesis
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
%%%%%%%%%%%%%%%%%%
if ~(strcmp(opt.noise,'homo')||strcmp(opt.noise,'hetero'))
    error('option ''opt.noise'' is wrongly specified!!')
end

% l
if numel(a.l) == 1
    a.l = [zeros(1,indx) a.l*ones(1,P-indx)];
else
    disp(['Veuillez donner ' num2str(P-indx) ' valeurs pour a.l']);
end

if numel(b.l) == 1
    b.l = [zeros(1,indx) b.l*ones(1,P-indx)];
else
    disp(['Veuillez donner ' num2str(P-indx) ' valeurs pour b.l']);
end

% beta2
if numel(a.beta2)==1
    a.beta2=a.beta2*ones(P,1);
end
if numel(b.beta2)==1
    b.beta2=b.beta2*ones(P,1);
end

% gamma2
if numel(a.gamma2) == K
    a.gamma2 = repmat(a.gamma2(:),[1 P-indx]);
    a.gamma2 = [zeros(K,indx) a.gamma2];
elseif numel(a.gamma2)==1
    a.gamma2= [zeros(K,indx) a.gamma2*ones(K,P-indx)];
elseif numel(a.gamma2)==(P-indx)
    a.gamma2= [zeros(K,indx) repmat(a.gamma2(:)',[K 1])];
elseif numel(a.gamma2)==K*(P-indx)
    a.gamma2=[zeros(K,indx) a.gamma2];
else
    disp(['Veuillez donner ' num2str(K) 'x' num2str(P-indx) ' valeurs pour a.gamma2']);
end
if numel(b.gamma2) == K
    b.gamma2 = repmat(b.gamma2(:),[1 P-indx]);
    b.gamma2 = [zeros(K,indx) b.gamma2];
elseif numel(b.gamma2) == 1
    b.gamma2 = [zeros(K,indx) b.gamma2*ones(K,P-indx)];
elseif numel(b.gamma2)==(P-indx)
    b.gamma2= [zeros(K,indx) repmat(b.gamma2(:)',[K 1])];
elseif numel(b.gamma2)==K*(P-indx)
    b.gamma2 = [zeros(K,indx) b.gamma2];
else
    disp(['Veuillez donner ' num2str(K) 'x' num2str(P-indx) ' valeurs pour b.gamma2']);
end

% sig2
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
%%%%%%%%%%%%%%%%%%
l = ones(Nrun,P); %no sparsity for the initialization
q = ones(Nrun,K,P);
q(1,:,:) = ones(K,P);
lq = K*ones(1,P); % initial number of non-nul binaries
beta2 = ones(Nrun,P); % default value is one

% can be initialized
Lambda = zeros(Nrun,M,K,P);
sig2 = zeros(Nrun,M);
gamma2 = zeros(Nrun,K,P);
Sc = zeros(Nrun,K,K,P);
s=ones(Nrun,K,P);

if nargin > 7
    if numel(Ini.Lambda)~=M*K*(P-indx); disp(['Veuillez donner ' num2str(M) 'x' num2str(K) 'x' num2str(P-indx) ' valeurs pour Ini.Lambda']); end;
    Lambda(1,:,:,1+indx:end) = Ini.Lambda; %(NrunxMxKxP)
    if numel(Ini.sig2)~=M; disp(['Veuillez donner ' num2str(M) ' valeurs pour Ini.sig2']); end;
    sig2(1,:) = Ini.sig2; %(NrunxM)
    if (numel(Ini.gamma2)~=K*(P-indx) && numel(Ini.gamma2)~=1); disp(['Veuillez donner ' num2str(K) 'x'  num2str(P-indx) ' valeurs pour Ini.gamma2']); end;
    gamma2(1,:,1+indx:end) = Ini.gamma2*ones(K,1); %(NrunxKxP)
    if isfield(Ini,'Sc')
        Sc(1,:,:,1+indx:end) = Ini.Sc;
    else
        for p=1+indx:P
            Sc(1,:,:,p) = cwishrnd(diag(gamma2(1,:,p)) , Isnap(p));
        end
    end
else
    Lambda(1,:,:,1+indx:end) = (randn(M,K,P-indx) + 1i*randn(M,K,P-indx))/sqrt(2*K_est); %(NrunxMxKxP)
    if opt.ref==1
        sig2(1,:)=real(diag(Sy(:,:,1)));
    else
        sig2(1,:) = b.sig2(:)'./(a.sig2(:)-1)';
    end
    gamma2(1,:,1+indx:end) = b.gamma2(1+indx:end)./(a.gamma2(1+indx:end)-1)*ones(K,1);
    for p=1+indx:P
        Sc(1,:,:,p) =  cwishrnd(diag(gamma2(1,:,p)) , Isnap(p));
    end
end

if isfield(opt,'gamma2')
    if opt.gamma2 > 0
        gamma2 = opt.gamma2*ones(Nrun,K,P);
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

regul=ones(1,Nrun);
%regul(1:round(Nrun/2))=logspace(log10(M*Isnap),log10(1),round(Nrun/2));
%regul(1:round(Nrun/2))=linspace(M*Isnap,1,round(Nrun/2));
% décroissance de la température en escalier
%Nmarches = 10;
%regul(1:round(Nrun/2))= ones(round(Nrun/2/Nmarches),1)*linspace(M*Isnap*10,1,Nmarches);



for i = 2:Nrun
    for p=1+indx:P
        
        Dq = diag(q(i-1,:,p));
        D=diag((beta2(i-1,p).*sig2(i-1,:)).^-1);
        
        % sample in [Sc|rest]
        % ===================
        LambdaD = squeeze2(Lambda(i-1,:,:,p))'*D;
        Cov_inv = Dq*LambdaD*squeeze2(Lambda(i-1,:,:,p))*Dq + diag(1./(gamma2(i-1,:,p)));
        Cov_inv = (Cov_inv + Cov_inv')/2;
        Cov_inv_sqrt = chol(Cov_inv);
        Cov_sqrt = Cov_inv_sqrt\eye(K);
        Cov = Cov_sqrt*Cov_sqrt';
        
        Wc = Cov_sqrt*cwishrnd([1 K],Isnap(p))*Cov_sqrt';
        temp = Cov*Dq*LambdaD*Sy(:,:,p)*LambdaD'*Dq*Cov' + Wc;
        Sc(i,:,:,p) = (temp + temp')/2;
        
        % sample in [Lambda|rest]
        % =======================
        Omega_inv = kron(Dq*conj(squeeze(Sc(i,:,:,p)))*Dq,sparse(D)) + lq(p)*speye(MK); %+K*speye(MK)
        Omega_inv = (Omega_inv + Omega_inv')/2;
        Omega_inv_sqrt = chol(Omega_inv);
        Omega_sqrt = inv(Omega_inv_sqrt);
        Omega = Omega_sqrt*Omega_sqrt';
        Syc = Sy(:,:,p)*D*squeeze(Lambda(i-1,:,:,p))*Dq*Cov;
        temp = D*Syc*Dq;
        temp = Omega*temp(:) + Omega_sqrt*(randn(MK,1) + 1i*randn(MK,1))/sqrt(2);
        Lambda(i,:,:,p) = reshape(temp,M,K);
        
        
        % Scaling of Sc and L
        %====================
        %{
        for k=1:K
            a_r = real(2/gamma2(i-1,k,p)*Sc(i,k,k,p));
            b_r = real(2*lq(p)* conj(Lambda(i,:,k,p))*Lambda(i,:,k,p).');
            p_r = Isnap(p) - M -0.5;
           
            r = gigrnd( p_r, a_r, b_r,1);
            
            s(i,k,p) = sqrt(r)*exp(1i*2*pi*rand); %(-1)^(round(rand)) * sqrt(r);
            Lambda(i,:,k,p) = Lambda(i,:,k,p)./s(i,k,p);
        end
        
        Lambda(i,:,k,p) = Lambda(i,:,k,p)./s(i,k,p);
        Sc(i,:,:,p) = squeeze(Sc(i,:,:,p)).*(s(i,:,p)'*s(i,:,p));
        
        %update Cov and Wc
        LambdaD = squeeze2(Lambda(i,:,:,p))'*D;
        Cov_inv = Dq*LambdaD*squeeze2(Lambda(i,:,:,p))*Dq + diag(1./(gamma2(i-1,:,p)));
        Cov_inv = (Cov_inv + Cov_inv')/2;
        Cov_inv_sqrt = chol(Cov_inv);
        Cov_sqrt = Cov_inv_sqrt\eye(K);
        Cov = Cov_sqrt*Cov_sqrt';
        Wc = Cov_sqrt*cwishrnd([1 K],Isnap(p))*Cov_sqrt';
        %}
        
        % sample in [gamma2|rest]
        % =======================
        if opt.gamma2 == 0
            %{
            a_gamma2 = a.gamma2(p) + K*Isnap(p);
            b_gamma2 = b.gamma2(p) + trace(squeeze(Sc(i,:,:,p)));
            gamma2(i,:,p) = ones(K,1)./rand_gamma(1,1,1/b_gamma2,a_gamma2);
            %}
                       
            a_gamma2 = a.gamma2(p) + Isnap(p);
            b_gamma2 = b.gamma2(p) + real(diag(squeeze(Sc(i,:,:,p))));
            for k=1:K
                gamma2(i,k,p) = 1./rand_gamma(1,1,1/b_gamma2(k),a_gamma2);
            end            
        end
        
        % sample in [l|rest]
        % ==================
        D=diag((beta2(i-1,p).*sig2(i-1,:)).^-1);
        a_l = a.l(p) + lq(p);
        b_l = b.l(p) + (K-lq(p));
        l(i,p) = betarnd(a_l,b_l);
        
        
        % sample in [q|rest] with or without marginalization on C
        % =======================================================
        for k=1:K           
            if strcmp(opt.marg,'off') %pas de marginalisation
                Z = eye(M) - squeeze(Lambda(i,:,:,p))*Dq*Cov*Dq*squeeze(Lambda(i,:,:,p))'*D;
                T = Z*Sy(:,:,p)*Z + squeeze(Lambda(i,:,:,p))*Dq*Wc*Dq*squeeze(Lambda(i,:,:,p))';
                %gini = real(trace(D*T)) + q(i-1,k,p).* log(1/l(i,p) -1) ;
                gini = real(trace(D*T)) + q(i-1,k,p).* log(1/l(i,p) -1) - vec(Lambda(i,:,:,p))'*vec(Lambda(i,:,:,p))*sum(q(i-1,:,p))+MK*log(sum(q(i-1,:,p)));
            elseif strcmp(opt.marg,'on') %avec marginalisation
                Bini = squeeze(Lambda(i,:,:,p))*Dq*diag(gamma2(i-1,:,p))*Dq*squeeze(Lambda(i,:,:,p))' + beta2(i-1,p)*diag(sig2(i-1,:));
                Bini=(Bini+Bini')/2;
                %gini = real(trace(inv(Bini)*Sy(:,:,p))) + Isnap(p)*log(real(det(Bini))) + q(i-1,k,p).* log(1/l(i,p) -1);
                gini = real(trace(inv(Bini)*Sy(:,:,p))) + Isnap(p)*log(real(det(Bini))) + q(i-1,k,p).* log(1/l(i,p) -1) - vec(Lambda(i,:,:,p))'*vec(Lambda(i,:,:,p))*sum(q(i-1,:,p))+MK*log(sum(q(i-1,:,p)));
            else 
                error('option.marg must be ''on'' or ''off''')                
            end
            dk = (-1)^q(i-1,k,p);
            Dqmod = Dq;
            Dqmod(k,k)= Dqmod(k,k) + dk;
            
            if strcmp(opt.marg,'off') %pas de marginalisation
                Z = eye(M) - squeeze(Lambda(i,:,:,p))*Dqmod*Cov*Dqmod*squeeze(Lambda(i,:,:,p))'*D;
                T = Z*Sy(:,:,p)*Z + squeeze(Lambda(i,:,:,p))*Dqmod*Wc*Dqmod*squeeze(Lambda(i,:,:,p))';
                %gmod =  real(trace(D*T)) + (q(i-1,k,p)+dk)* log(1/l(i,p) - 1);
                gmod =  real(trace(D*T)) + (q(i-1,k,p)+dk)* log(1/l(i,p) - 1) -  vec(Lambda(i,:,:,p))'*vec(Lambda(i,:,:,p))*sum(diag(Dqmod))+MK*log(sum(diag(Dqmod)));
            elseif strcmp(opt.marg,'on')  %avec marginalisation
                Bmod = Bini + dk * Lambda(i,:,k,p).'*gamma2(i-1,k,p)*conj(Lambda(i,:,k,p));
                Bmod=(Bmod+Bmod')/2;
                %gmod = real(trace(inv(Bmod)*Sy(:,:,p))) + Isnap(p)*log(real(det(Bmod))) + (q(i-1,k,p)+dk).* log(1/l(i,p) -1); % scalar
                gmod = real(trace(inv(Bmod)*Sy(:,:,p))) + Isnap(p)*log(real(det(Bmod))) + (q(i-1,k,p)+dk).* log(1/l(i,p) -1) -  vec(Lambda(i,:,:,p))'*vec(Lambda(i,:,:,p))*sum(diag(Dqmod))+MK*log(sum(diag(Dqmod)));
            else 
                error('option.marg must be ''on'' or ''off''')
            end
            g = real(gini-gmod);
            %g = real(gini(k)-gmod);
            
            temp = rand;
            %eq. 4.81 C. FAURE thesis
            
            %if temp < 1/(1 + exp(-g))
            if log(1/(temp) - 1)> -g
                q(i,k,p) =  q(i-1,k,p) + dk; % le binaire change d'état par rapport au tirage de l'itération précédente
            else
                q(i,k,p) =  q(i-1,k,p);
            end
            %q(i,k,p)=1;
            %Dq = diag([q(i,1:k,p) q(i-1,(k+1):K,p)]);
        end
        indq{p} = find(q(i,:,p)==1); %indices in 1:K for which q is 1
        lq(p) = length(indq{p});
        
        %{
        if lq(p)==0
            q(i,1,p)=1;
            lq(p)=1;
            indq{p}=1;
        end
        %}
        Dq = diag(q(i,:,p));
        
        
        % sample in [beta2|rest]
        % ======================
        Z = eye(M) - squeeze(Lambda(i,:,:,p))*Dq*Cov*Dq*squeeze(Lambda(i,:,:,p))'*D;
        for k = 1:M
            Er(k,p) = real(Z(k,:)*Sy(:,:,p)*Z(k,:)' + squeeze(Lambda(i,k,:,p)).'*Dq*Wc*Dq*conj(squeeze(Lambda(i,k,:,p))));
        end
        b_beta2 = b.beta2(p) + sig2(i-1,:).^-1*Er(:,p);
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        beta2(i,p) = 1;%1/rand_gamma(1,1,1/b_beta2,a_beta2(p));
        
    end
    
    % sample in [sig2|rest]
    % =======================
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
%close(h)

end


