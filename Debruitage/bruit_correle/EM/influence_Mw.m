%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all
addpath(genpath('../../..'))
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

freq = 3800;
%Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),50));
Nsrc =25;
rho=0;
SNR=10;

option.max = 3000; %max number of iteration
option.rerr=1e-6;


for i=1:length(Mw)
	i    
    if Mw(i)<93
        k=Mw(i)-1;
    else
        k=92;
    end
    
    %%% Generate data
    [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw(i));
    %Sy(:,:)=Sp+((Sn));
    %Rang(i,j)=rank(Sp,0.005*max(eig(Sp)));
    Rang2(i)=rank(Sp);
    
    CSM = Sy(:,:);
    d_ref(:,i)=real(diag(Sy-Sn)); 
    

	 
	%initialisation nulle
    Ini.Syc = 1e-16*ones(93,1);
    Ini.L=1e-16*ones(93,k);
    
    [L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,93,k,option,Ini);
    flag.count    
    
    d_EM(:,i)=real(diag(Sx));

    norm_k(i) = norm( d_EM(:,i)- d_ref(:,i)) /  norm( d_ref(:,i));
    


    
    norm_it_rang(1,i)=norm(real(diag(Ini.L*Ini.L'))-d_ref(:,i))/norm(d_ref(:,i));
    for j=1:flag.count
        norm_it_rang(j+1,i)=norm(real(d1all(:,j))-d_ref(:,i))/norm(d_ref(:,i));
        %norm_it_rang2(j+1,i)=norm(real(d2all(:,j))-d_ref(:,i))/norm(d_ref(:,i));
    end
   
end

%save('EM_Rang','norm_k','Rang2','d_EM','d_ref');
