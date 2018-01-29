
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =1:93;
Mw=10^4;
rho=0;
SNR =10;
M=93; %nb of microphones


j=1;

Nrun=1000;
opt.noise='hetero';
K_est=93;
alpha_save=zeros(M);

for i=1:length(Nsrc)
	disp(['nb sources : ' num2str(Nsrc(i))])
    
    K_est=Nsrc(i)+5;
    if K_est>93
        K_est=93;
    end
    K_est
	
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw);
	d_ref(:,i,j)=real(diag(Sp));
      
	%%% initialisation
	
	sign = 1*(1 + .1*rand(M,1));	
	kappa = 0;      % important!!
	%alpha2_mean = exp(-kappa*(0:K_est-1)'/(K_est-1));
    alpha2_mean = linspace(max(real(eig(Sy))),min(real(eig(Sy))),M)./max(real(eig(Sy)));
	
	%alpha2
%	for k = 1:K_est
%		[a.alpha2(k),b.alpha2(k)] = Convert_InvGamma(alpha2_mean(k),100*alpha2_mean(k));   % hyper-hyper-paramètres sur alpha2
%	end

	%alpha
	for k = 1:K_est
    	a.alpha(k) = 1/alpha2_mean(k);      % hyper-hyper-paramètres sur alpha
	end	
	
	%beta
	[a.beta2,b.beta2] = Convert_InvGamma(mean(sign.^2),10*mean(sign.^2));     % hyper-hyper-paramètres sur beta2
	
	%gamma
	gamma_mean = (real(trace(Sy))/M )/mean(alpha2_mean); %Alice : pas de bruit retiré
	[a.gamma2,b.gamma2] = Convert_InvGamma(gamma_mean,10*gamma_mean);   % hyper-hyper-paramètres sur gamma2
	% option.gamma2 = 100;

    a.alpha2=a.alpha;
    b.alpha2=1.1;
    %opt='hetero';
	[Sc,Lambda,alpha,beta2] = MCMC_AnaFac_Quad_Sparse3(Sy,K_est,a,b,Mw,Nrun,opt);
	alpha_save(1:K_est,i)=alpha(end,:);
	
	for jj=1:Nrun
		d_mcmc(:,jj,i)=real(diag(squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))'));
        err_mcmc(i,jj) = norm( d_ref(:,i) - d_mcmc(:,jj,i) ) /  norm( d_ref(:,i));
	end
	

% 	norm_it_rang(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref(:,i))/norm(d_ref(:,i));
%	for jj=1:flag.count
%		norm_relative(jj+1,i)=norm(real(d1all(:,jj))-d_ref(:,i))/norm(d_ref(:,i));
%		norm_sig2(jj+1,i)=flag.norm(jj);
%		loglik_all(jj+1,i)=loglik(jj);
%		diffloglik_all(jj+1,i)=diffloglik(jj);
%		%norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
%	end

end

save('MCMC_rang','err_mcmc','SNR','Nsrc','d_mcmc', 'd_ref','alpha_save');
%figure
%plot(Nsrc,10*log10(err_EM))



