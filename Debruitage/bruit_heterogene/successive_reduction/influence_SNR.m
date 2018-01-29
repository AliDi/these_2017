%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =20;
Mw=10^4;
rho=0;
SNR = -10:1:10;
%extra_SNR=0;

for i=1:length(SNR)
i
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal_hetero(freq, Nsrc, rho , SNR(i) , Mw,SNR(i)-10);

	CSM = Sy(:,:);

	d_ref(:,i)=real(diag(Sp));


	%%%--------------------------------------------------------------------------------------------
	%%% Alternating projections
	%%%--------------------------------------------------------------------------------------------
	%cvx (Hald)
	cvx_quiet('true')
	cvx_precision('high')
	[CSM_cvx d1 cvx_it]=CSMRecHald(CSM);     
	d_cvx(:,i)=diag(CSM_cvx);
	err_cvx(i)= norm( d_ref(:,i)-real(d_cvx(:,i)) ,2)/norm(d_ref(:,i),2); 


	 %linprog (Dougherty)
     [CSM_linprog ii] = recdiagd(double(CSM),500,60);
     d_linprog(:,i) = real(diag(CSM_linprog));
     err_linprog(i)= norm( d_ref(:,i)-d_linprog(:,i)) / norm(d_ref(:,i)) ;

	%Alternating projection
    [CSM_it, n, errec, D_AP] = recdiag(CSM,1,1000,1e-7,30); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
    d_it(:,i) = real(diag(CSM_it));
    err_it(i) = norm( d_ref(:,i)-d_it(:,i)) / norm(d_ref(:,i));

	%Version Jérôme (remove average value of the K greatest SV)
	%[Sn,L] =
%         SS_CSM_Fit(CSM,92);
	%d_jerome(:,i) = real(diag(L*L'));
	%d_jerome2(:,i) = real(diag(CSM-diag(Sn)));
	%err_jerome(i) = norm( d_ref(:,i)-real(d_jerome(:,i))) / norm(d_ref(:,i));

end

save('cvx_SNR','err_cvx','SNR','Nsrc','d_cvx', 'd_ref');
save('linprog_SNR','err_linprog','SNR','d_linprog','d_ref');
save('AP_it_SNR','err_it','SNR','Nsrc','d_it','d_ref');







