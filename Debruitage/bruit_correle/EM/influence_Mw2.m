%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),50));
rho=0;
SNR = 10;

option.max = 1000; %max number of iteration
option.rerr=1e-7;


for j=1:length(Mw)
    j
	for i=1:length(Nsrc)
	i
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j));
		Rang2(i,j)=rank(Sp);
		CSM = Sy(:,:);
		d_ref(:,i,j)=real(diag(Sy-Sn));
		
		
		if Mw(j)<93
       		 k=Mw(j)-1;
    	else
        	k=92;
    	end
		
		%initialisation nulle
		Ini.Syc = 1e-16*ones(93,1);
		Ini.L=1e-16*ones(93,k);
		
		    
		[L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,93,k,option,Ini);
		flag.count    
		
		d_EM(:,i,j)=real(diag(Sx));
		err_EM(i,j) = norm( d_ref(:,i,j) - d_EM(:,i,j) ) /  norm( d_ref(:,i,j));

	end
end
save('EM_Mw','err_EM','Mw','Nsrc','d_EM', 'd_ref');




