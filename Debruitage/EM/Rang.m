%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all
addpath(genpath('..'))
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

freq = 3800;
Mw=9000;
Nsrc =1:2:150;
rho=0;
SNR=10;

option.max = 2000; %max number of iteration
option.rerr=5e-5;

k=92;

for i=1:length(Nsrc)
	i
    %%% Generate data
    [Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw);
    Sy(:,:)=Sp+diag(diag((Sn)));
    %Rang(i,j)=rank(Sp,0.005*max(eig(Sp)));
    Rang2(i)=rank(Sp);
    
    CSM = Sy(:,:);
    d_ref(:,i)=real(diag(Sp));
    
    %[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy,k); %initialisation SI
	 
	%initialisation nulle
    Ini.Syc = 1e-16*ones(93,1);
    Ini.L=1e-16*ones(93,k);
    
    [L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,Mw,k,option,Ini);
    flag.count    
    
    d_EM(:,i)=real(diag(Sx));
    norm_k(i) = norm( d_EM(:,i)- d_ref(:,i)) /  norm( d_ref(:,i));
   
end

save('EM_Rang','norm_k','Rang2','d_EM','d_ref');
