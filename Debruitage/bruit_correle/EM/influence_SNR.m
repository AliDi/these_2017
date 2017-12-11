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
Mw=10^4;
rho=0;
SNR = -10:10;

option.max = 1000; %max number of iteration
option.rerr=1e-7;

for j=1:length(Mw)
    j
	for i=1:length(SNR)
	i
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw);
		
		CSM = Sy(:,:);
		d_ref(:,i,j)=real(diag(Sp));
				
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
		iter(i)=flag.count;  
		
		d_EM(:,i,j)=real(diag(Sx));
		err_EM(i,j) = norm( d_ref(:,i,j) - d_EM(:,i,j) ) /  norm( d_ref(:,i,j));
		
%		 norm_it_rang(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref(:,i))/norm(d_ref(:,i));
%		for jj=1:flag.count
%		    norm_it_rang(jj+1,i)=norm(real(d1all(:,jj))-d_ref(:,i))/norm(d_ref(:,i));
%		    %norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
%		end

	end
end
save('EM_SNR','err_EM','SNR','Nsrc','d_EM', 'd_ref');




