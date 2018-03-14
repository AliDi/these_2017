
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =1:93;
Mw=10^4;

rho=0;
SNR = 10;
M=93; %nb of microphones


K_est=M;
Nrun=1000;
opt.noise='hetero';

d_mcmc=zeros(M,Nrun,length(Nsrc));
err_mcmc=zeros(length(Nsrc),Nrun);
alpha_save=zeros(M,length(Nsrc));

parfor i=1:length(Nsrc)
	disp(['Nsrc : ' num2str(Nsrc(i))])
        
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw , SNR );

	d_ref=real(diag(Sp));
    save_d_ref(:,i)=d_ref;
			
    
    %%% initialisation	
    
    a=struct();
    b=struct();
    Ini=struct();
    
	noise = real(mean(diag(Sy)));%1*(1 + .1*rand(M,1));	
	kappa = 1;      % important!!
	%alpha2_mean = exp(-kappa*(0:K_est-1)'/(K_est-1));
    alpha2_mean = abs(real(sort(eig(Sy)/max(eig(Sy)),'descend')));
	   
    %%% hyper-parametres
	%alpha2
%	for k = 1:K_est
%		[a.alpha2(k),b.alpha2(k)] = Convert_InvGamma(alpha2_mean(k),100*alpha2_mean(k));   % hyper-hyper-paramètres sur alpha2
%	end

	%alpha
	for k = 1:K_est
    	a.alpha(k) = 1/alpha2_mean(k);      % hyper-hyper-paramètres sur alpha
	end	
	
	%beta
	[a.beta2,b.beta2] = Convert_InvGamma(mean(noise.^2),10*mean(noise.^2));     % hyper-hyper-paramètres sur beta2
	
	%gamma
	gamma_mean = (real(trace(Sy))/M )/mean(alpha2_mean); %Alice : pas de bruit retiré
	[a.gamma2,b.gamma2] = Convert_InvGamma(gamma_mean,10*gamma_mean);   % hyper-hyper-paramètres sur gamma2
	% option.gamma2 = 100;
    
    %initialisation
    Ini.alpha(1,:) = 1./a.alpha(:)';
    Ini.beta2(1,:) = b.beta2/a.beta2*ones(M,1);
    Ini.gamma2(1) = b.gamma2/a.gamma2;
    Ini.Lambda(1,:,:) = (randn(M,K_est) + 1i*randn(M,K_est))/sqrt(2);
    
    %opt.gamma2=1.1;
    

    
	[Sc,Lambda,alpha,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3(Sy,K_est,a,b,Mw,Nrun,opt,Ini);
    
    tmp_alpha_save=zeros(M,1);
    tmp_alpha_save(1:K_est)=alpha(end,:);
	alpha_save(:,i)=tmp_alpha_save;
	
    tmp_d_mcmc=zeros(M,Nrun);
    tmp_err_mcmc=zeros(Nrun,1);
	for jj=1:Nrun
		tmp_d_mcmc(:,jj)=real(diag(squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))'));
        tmp_err_mcmc(jj,1) = norm( d_ref - tmp_d_mcmc(:,jj)) /  norm( d_ref);        
    end
    d_mcmc(:,:,i)=tmp_d_mcmc;
    err_mcmc(i,:)=tmp_err_mcmc;
	
end

save('MCMC_Nsrc_alphaVP','err_mcmc','Mw','Nsrc','d_mcmc', 'save_d_ref','alpha_save');
%figure
%plot(Nsrc,10*log10(err_EM))



