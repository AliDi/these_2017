%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

freq = 3000;
Nsrc = 100;
rho=0;
SNR=-10:2:10;

for snr=1:length(SNR)

	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(snr));

	CSM = Sy;
	d_ref = real(diag(Sp));

	%%%--------------------------------------------------------------------------------------------
	%%% SLDR solved with proximal gradient
	%%%--------------------------------------------------------------------------------------------
	lambda= 0.1;
	[A E Nit out A_all] = proximal_gradient_rpca(CSM , lambda, 70,1e-10,-1,-1,-1,-1);
	
	d_pg = real( diag(A_all(:,:,end)));	
	err_pg(snr)= norm( diag(d_ref - d_pg) ) / norm(d_ref);

	%%%--------------------------------------------------------------------------------------------
	%%% Alternating projections
	%%%--------------------------------------------------------------------------------------------
	%cvx (Hald)
	cvx_quiet('true')
	[CSM_cvx d1 cvx_it]=CSMRecHald(CSM);
	
	d_cvx = real(diag(CSM_cvx));
	err_cvx(snr)= norm( diag(d_ref - d_cvx) ) / norm(d_ref);

	%linprog (Dougherty)
	%[CSM_linprog ii x1 D_linprog(:,:)] = recdiagd(double(CSM),1000,45);

	%Alternating projection
	%[CSM_it, n, errec, D_AP] = recdiag(CSM,1,70,1e-8,60); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
end






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

