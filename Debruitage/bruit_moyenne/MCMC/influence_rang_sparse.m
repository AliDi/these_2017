
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =20%[20 80];%1:93;
Mw=10^4;
rho=0;
SNR = 10;
M=93; %nb of microphones

j=1;
for i=1:length(Nsrc)
disp(['nb sources : ' num2str(Nsrc(i))])
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw);

	d_ref(:,i,j)=real(diag(Sp));
			
    k=Nsrc(i)+round(0.5*Nsrc(i));
    if k>M
        k=M;
    end
	%initialisation
	
	%Ini.Lambda=(randn(M,k) + 1i*randn(M,k))/sqrt(2);
	sign = 1*(1 + .1*rand(M,1));	
	[a.alpha2,b.alpha2] = Convert_InvGamma(1,100);   % hyper-hyper-paramètres sur alpha2
	[a.beta2,b.beta2] = Convert_InvGamma(sign(1)^2,10*sign(1)^2);     % hyper-hyper-paramètres sur beta2
	[a.gamma2,b.gamma2] = Convert_InvGamma(1,10);   % hyper-hyper-paramètres sur alpha2
    
   % Ini.alpha2= 1/rand_gamma(1,1,1/b.alpha2,a.alpha2);
    %Ini.beta2= 1/rand_gamma(1,1,1/b.beta2,a.beta2);
    		
	Nrun=100;
	opt='hetero';
	
	
	[Sc,Lambda,alpha2,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse(Sy,k,a,b,Mw,Nrun,opt);
	%[Sc,Lambda,alpha2,beta2] = MCMC_AnaFac_Quad(Sy,k,a,b,Mw,Nrun,opt);
	
	for jj=1:Nrun
		d_mcmc(:,jj,i)=real(diag(squeeze(Lambda(jj,:,:))*squeeze(Sc(jj,:,:))*squeeze(Lambda(jj,:,:))'));
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

%save('test_MCMC_Nsrc','err_EM','SNR','Nsrc','d_EM', 'd_ref');
%figure
%plot(Nsrc,10*log10(err_EM))



