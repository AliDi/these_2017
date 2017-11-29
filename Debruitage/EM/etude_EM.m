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
Nsrc =25;
%Nsrc=ceil(sqrt(Nsrc)).*floor(sqrt(Nsrc));
%Nsrc=unique(Nsrc);
rho=0;
SNR=10;

option.max = 2000; %max number of iteration
option.rerr=1e-7;

[Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw);
Sy(:,:)=Sp+diag(diag((Sn)));
d_ref = real(diag(Sp));


for k=1:92
k
    %[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy,k); %initialisation
	  
    Ini.Syc = 1e-16*ones(93,1);
    Ini.L=1e-16*ones(93,k);
    
    [L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,Mw,k,option,Ini);
    flag.count    
    
    norm_k(k) = norm( real(diag(Sx)) - d_ref ) /  norm( d_ref);
    
    norm_it_k(1,k)=norm(real(diag(Ini.L*Ini.L'))-d_ref)/norm(d_ref);
    for i=1:flag.count
        norm_it_k(i+1,k)=norm(real(d1all(:,i))-d_ref)/norm(d_ref);
        %norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
    end


    
    
end


