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
SNR = 10;

option.max = 200; %max number of iteration
option.rerr=0.5e-3;

K=1:92;

%%% Generate data
[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw);

CSM = Sy(:,:);
d_ref=real(diag(Sp));

j=1;
for i=1:length(K)
i	
	%initialisation nulle
	%Ini.Syc = 1e-16*ones(93,1);
	%Ini.L=1e-16*ones(93,K(i));
	
	[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy,K(i));
	
	    
	[L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,93,K(i),option,Ini);
	flag.count  
	iter(i)=flag.count;  
	
	d_EM(:,i,j)=real(diag(Sx));
	err_EM(i,j) = norm( d_ref - d_EM(:,i,j) ) /  norm( d_ref);
	
	 norm_it_rang(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref)/norm(d_ref);
	for jj=1:flag.count
	    norm_it_rang(jj+1,i)=norm(real(d1all(:,jj))-d_ref)/norm(d_ref);
	    %norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
	end
			%a(i)=sum(sig2)
			%size(sig2)
end

save('test_EM_K','err_EM','K','Nsrc','d_EM', 'd_ref');
plot(K,10*log10(err_EM))



