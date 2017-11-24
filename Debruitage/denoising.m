%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | 
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

%%%--------------------------------------------------------------------------------------------
%%% Load CSM signal
%%%--------------------------------------------------------------------------------------------
Spp_filename = 'generate_spectra/Sp_SNR0dB_Mw9300_100src.h5';

Sp = h5read(Spp_filename , '/signal');
freqs = h5read(Spp_filename , '/frequencies');
Nfreq=length(freqs);
Sy = h5read(Spp_filename, '/signal_plus_noise');
Sn = h5read(Spp_filename, '/noise');

Nmic=size(Sp,1);
Mw=9300;

f=3;

d_ref = real(diag(Sy(:,:,f)-Sn(:,:,f)));
CSM = Sy(:,:,f);

%%

%%%--------------------------------------------------------------------------------------------
%%% SLDR solved with proximal gradient
%%%--------------------------------------------------------------------------------------------
lambda= 0.1;

[A E Nit out A_all] = proximal_gradient_rpca(CSM , lambda, 70,1e-10,-1,-1,-1,-1);

for i=1:Nit
	norm_pg(i) = norm((diag(A_all(:,:,i))-d_ref)./d_ref);
end

%%DRfinez
%CSM_DRfinez=recdiagf(double(CSM));%hypothese de coherence à verifier



%%%--------------------------------------------------------------------------------------------
%%% Alternating projections
%%%--------------------------------------------------------------------------------------------

%cvx (Hald)
precision=logspace(log10(1-1e-10),log10(1e-13),100);
cvx_it_former=0;
D_cvx=[];
for i=1:100
	cvx_quiet('true')
	cvx_precision([precision(i) precision(i) precision(i)]);
	[CSM_cvx d1 cvx_it(i)]=CSMRecHald(CSM);
	if cvx_it_former<cvx_it(i)
		D_cvx=[D_cvx diag(CSM_cvx)];
	end
	cvx_it_former=cvx_it(i);
end

cvx_it=unique(cvx_it);

for i=1:length(cvx_it)
    norm_cvx(i)=norm((D_cvx(:,i)-d_ref)./d_ref);
end


%linprog (Dougherty)
D_linprog=[];
[CSM_linprog ii x1 D_linprog(:,:)] = recdiagd(double(CSM),1000,45);
for i=1:ii
    norm_linprog(i)=norm((D_linprog(:,i)-d_ref)./d_ref);
end

%Alternating projection
[CSM_it, n, errec, D_AP] = recdiag(CSM,1,70,1e-8,60); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
for i=1:n
    norm_AP(i)=norm((D_AP(:,i)-d_ref)./d_ref);
end


%%%--------------------------------------------------------------------------------------------
%%% Keep only K greatest eigenvalues
%%%--------------------------------------------------------------------------------------------
s= min(real(eig(CSM)))/(1-sqrt(Nmic/Mw))^2; %estimation du bruit
[l cdf]=distrib_vp_noise(Nmic,Mw,s);
figure
plot(cdf,10*log10(l));
hold on
plot(10*log10(sort(real(eig(Sy(:,:,f))))));
K=input('Rang supposé ? ');

[noise, L_signal] = SS_CSM_Fit(CSM,K);
CSM_KVP = L_signal*L_signal';

%%%--------------------------------------------------------------------------------------------
%%% EM
%%%--------------------------------------------------------------------------------------------


%EM Jérôme
	%Initialisation par conservation des K premières VP
[Ini.Syc, Ini.L] = SS_CSM_Fit(Sy(:,:,f),K); %initialisation
    %Initialisation par des valeurs aléatoires
v=var(diag(CSM));
m=mean(diag(CSM));
Ini.Syc=v*randn(93,1)+m;
[u s] = eig(CSM-diag(Ini.Syc));
Ini.L=u(:,1:K)*sqrt(s(1:K,1:K));


option.max = 70; %max number of iteration
option.rerr=1e-7;

[L,sig2,beta2,flag,Sx, D_EM, D2] = EM_CSM_Fit(Sy(:,:,f),Mw,K,option,Ini,Sp(:,:,f));

for i=1:flag.count
    norm_EM(i)=norm(real(D_EM(:,i))-d_ref)./norm(d_ref);
    norm_EM2(i)=norm(real(D2(:,i))-d_ref)./norm(d_ref);
end

%display error
figure
hold on
plot(cvx_it,10*log10(norm_cvx))
plot(10*log10(norm_linprog))
plot(10*log10(norm_AP))
plot(10*log10(norm_pg))
plot(10*log10(norm_EM))
legend('CVX','Linprog','AP','PG','EM')

%diaplay denoised diagonal
figure
hold on
plot(diag(CSM));
plot(d_ref,'o-')
plot(D_cvx(:,end))
plot(D_linprog(:,end))
plot(D_AP(:,end))
plot((diag(A_all(:,:,end))))
plot(D_EM(:,end))
plot(diag(CSM_KVP))
legend('Noisy diagonal','Objective','CVX','Linprog','AP','PG','EM','K VP')

figure
hold on
plot(d_ref,diag(CSM),'o');
plot(d_ref,d_ref,'-')
plot(d_ref,D_cvx(:,end),'o')
plot(d_ref,D_linprog(:,end),'o')
plot(d_ref,D_AP(:,end),'o')
plot(d_ref,(diag(A_all(:,:,end))),'o')
plot(d_ref,D_EM(:,end),'o')
plot(d_ref,diag(CSM_KVP),'o')
legend('Noisy diagonal','Objective','CVX','Linprog','AP','PG','EM','K VP')

figure
hold on
plot((D_cvx(:,end)-d_ref).^2)
plot((D_linprog(:,end)-d_ref).^2)
plot((D_AP(:,end)-d_ref).^2)
plot((diag(A_all(:,:,end))-d_ref).^2)
plot((D_EM(:,end)-d_ref).^2)
plot((diag(CSM_KVP)-d_ref).^2)
legend('CVX','Linprog','AP','PG','EM','K VP')


%figure; hold on
%plot(real(diag(CSM)))
%%title('raw csm')

%plot(real(diag(CSM_DRfinez)))
%%title('csm finez')

%plot(real(diag(CSM_SLRD)))
%%title('csm sldr')

%plot(real(diag(CSM_cvx)))
%%title('csm cvx')

%plot(real(diag(CSM_linprog)))
%%title('csm linprog')


%plot(real(diag(CSM_it)))
%%title('csm iterative suppression')

%legend('raw','finez','slrd','cvx','linprog','AP');

