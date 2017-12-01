%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all
addpath(genpath('../..'))
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

freq = 3800;
Mw=9000;
Nsrc =5;
rho=0;
SNR=-20:2:20;

option.max = 1000; %max number of iteration
option.rerr=1e-8;

k=92;

for i=1:length(SNR)
	i
    %%% Generate data
    [Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw);
    Sy(:,:)=Sp+diag(diag((Sn)));
    %Rang(i,j)=rank(Sp,0.005*max(eig(Sp)));
    Rang2(i)=rank(Sp);
    
    CSM = Sy(:,:);
    d_ref(:,i)=real(diag(Sp));
    
    %[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy,k); %initialisation SI
	 
	%initialisation nulle
    Ini.Syc = 1e-16*ones(93,1);
    Ini.L=1e-16*ones(93,k);
    
    %[Ini.Syc, Ini.L] = SS_CSM_Fit(CSM,92); %initialisation
    
    
    [L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,93,k,option,Ini);
    flag.count    
    
    d_EM(:,i)=real(diag(Sx));
    norm_k(i) = norm( d_ref(:,i) - d_EM(:,i) ) /  norm( d_ref(:,i));
    %norm_k2(i)
    
    norm_it_snr(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref(:,i))/norm(d_ref(:,i));
    for j=1:flag.count
        norm_it_snr(j+1,i)=norm(real(d1all(:,j))-d_ref(:,i))/norm(d_ref(:,i));
        %norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)/norm(d_ref);
    end
   
end

%save('EM_SNR_5src','norm_k','SNR','d_EM','d_ref');
