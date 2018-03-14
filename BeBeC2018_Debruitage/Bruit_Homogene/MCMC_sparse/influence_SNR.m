
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =20;
Mw=10^4;
rho=0;
SNR = 0%-10:10;
M=93; %nb of microphones

j=1;
K_est=93;
Nrun=1000;
opt.noise='hetero';

for i=1:length(SNR)
	disp(['SNR : ' num2str(SNR(i))])
	
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw , SNR(i) );
	d_ref(:,i,j)=real(diag(Sp));
			
    
	%%% initialisation
	
	sign = 1*(1 + .1*rand(M,1));	
	kappa = 0;      % important!!
	%alpha2_mean = exp(-kappa*(0:K_est-1)'/(K_est-1));
    alpha2_mean = linspace(max(real(eig(Sy))),min(real(eig(Sy))),M)./max(real(eig(Sy)));%real(sort(eig(Sy)/max(eig(Sy)),'descend'));
	
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


	[Sc,Lambda,alpha,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3(Sy,K_est,a,b,Mw,Nrun,opt);
	
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

%save('MCMC_SNR','err_mcmc','SNR','Nsrc','d_mcmc', 'd_ref');
%figure
%plot(Nsrc,10*log10(err_EM))



