%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all
addpath(genpath('.'))
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

freq = 3000;
Mw=10000;
Nsrc =25;
%Nsrc=ceil(sqrt(Nsrc)).*floor(sqrt(Nsrc));
%Nsrc=unique(Nsrc);
rho=0;
SNR=-10;

option.max = 100; %max number of iteration
option.rerr=1e-6;

[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw);
d_ref = real(diag(Sp));


for k=1:92
    [Ini.Syc, Ini.L] = SS_CSM_Fit(Sy,k); %initialisation
%     v=var(diag(Sy));
%     m=mean(diag(Sy));
%     Ini.Syc=v*randn(93,1)+m;
%     [u s] = eig(Sy-diag(Ini.Syc));
%     Ini.L=u(:,1:k)*sqrt(s(1:k,1:k));
    
    [L,sig2,beta2,flag,Sx, d1all, d2all] = EM_CSM_Fit(Sy,Mw,k,option,Ini,Sp);
    
    %norm_EM1(k) = norm( real(diag(L*L')) - d_ref ,2) /  norm( d_ref,2);
    %norm_EM2(k) = norm( real(diag(Sx)) - d_ref ,2) /  norm( d_ref,2);
    
    for i=1:flag.count
        norm_EM(i,k)=norm(real(d1all(:,i))-d_ref)./norm(d_ref);
        norm_EM2(i,k)=norm(real(d2all(:,i))-d_ref)./norm(d_ref);
    end


    
    
end
%%

%%%--------------------------------------------------------------------------------------------
%%% Keep only K greatest eigenvalues
%%%--------------------------------------------------------------------------------------------
%s= min(real(eig(CSM)))/(1-sqrt(Nmic/Mw))^2; %estimation du bruit
%[l cdf]=distrib_vp_noise(Nmic,Mw,s);
%figure
%plot(cdf,10*log10(l));
%hold on
%plot(10*log10(sort(real(eig(Sy(:,:,f))))));
%K=input('Rang supposé ? ');

%[noise, L_signal] = SS_CSM_Fit(CSM,K);
%CSM_KVP = L_signal*L_signal';

%%%--------------------------------------------------------------------------------------------
%%% EM
%%%--------------------------------------------------------------------------------------------

%EM Jérôme
	%Initialisation par conservation des K premières VP
%[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy(:,:,f),K); %initialisation
%    %Initialisation par des valeurs aléatoires
%v=var(diag(CSM));
%m=mean(diag(CSM));
%Ini.Syc=v*randn(93,1)+m;
%[u s] = eig(CSM-diag(Ini.Syc));
%Ini.L=u(:,1:K)*sqrt(s(1:K,1:K));


%option.max = 70; %max number of iteration
%option.rerr=1e-7;

%[L,sig2,beta2,flag,Sx, D_EM] = EM_CSM_Fit(Sy(:,:,f),Mw,K,option,Ini,Sp(:,:,f));

%for i=1:flag.count
%    norm_EM(i)=norm((D_EM(:,i)-d_ref)./d_ref);
%end

