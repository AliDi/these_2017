%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =20;
Mw=round(logspace(log10(10),log10(50000),40));
rho=0;
SNR = 10;
extra_SNR=0;

option.max = 200; %max number of iteration
option.rerr=0.5e-8;
j=1;
for i=1:length(Mw)
i
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw(i),extra_SNR);
	
	CSM = Sy(:,:);
	d_ref(:,i,j)=real(diag(Sp));
			
	if Mw(i)<93
   		 k=Mw(i)-1;
	else
    	k=93;
	end
	
	%initialisation nulle
	Ini.Syc = 1e-16*ones(93,1);
	Ini.L=1e-16*ones(93,k);
	
	    
	[L,sig2,beta2,flag,Sx, d1all, d2all,loglik,diffloglik,sig2_all] = EM_CSM_Fit(Sy,Mw(i),k,option,Ini);
	flag.count  
	iter(i)=flag.count;  

	d_EM(:,i)=real(diag(Sx));
	err_EM(i) = norm( d_ref(:,i) - d_EM(:,i) ) /  norm( d_ref(:,i));

 	norm_it_rang(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref(:,i))/norm(d_ref(:,i));
	for jj=1:flag.count
		norm_relative(jj+1,i)=norm(real(d1all(:,jj))-d_ref(:,i))/norm(d_ref(:,i));
		norm_sig2(jj+1,i)=flag.norm(jj);
		loglik_all(jj+1,i)=loglik(jj);
		diffloglik_all(jj+1,i)=diffloglik(jj);
		sig2_all_all(jj+1,i)=sig2_all(jj);
		%norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
	end

end
save('EM_Mw','err_EM','Mw','Nsrc','d_EM', 'd_ref');
%figure
%semilogx(Mw,10*log10(err_EM))



