%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 3800;
Nsrc =1:1:150;
Mw=9000;
rho=0;
SNR = 10;

for i=1:length(Nsrc)
i
	%%% Generate data
	[Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw);
	Sy(:,:,i)=Sp+diag(diag((Sn)));
	%Rang(i,j)=rank(Sp,0.005*max(eig(Sp)));
	Rang2(i)=rank(Sp);

	CSM = Sy(:,:,i);
	d_ref(:,i)=real(diag(Sp));

	%%%--------------------------------------------------------------------------------------------
	%%% RPCA solved with proximal gradient
	%%%--------------------------------------------------------------------------------------------
	lambda= 0.5;
	[A E Nit out] = proximal_gradient_rpca(CSM , lambda, 300,1e-8,-1,-1,-1,-1);

	d_pg(:,i) = real( diag(A) );	
	err_pg(i)= norm( diag(d_ref(:,i) - d_pg(:,i))) / norm(d_ref(:,i));
end

save('RPCA_rang','err_pg','Rang2','d_pg','d_ref');

[c ia ic]=unique(Rang2);

figure
plot(c,10*log10(err_pg(ia)))
xlabel('Rank of $\boldmath{S_p}$')
ylabel('Relative error on diag($\boldmath{S_y}$) (dB)')
xlim([0 93])
plot_fig(gcf,11,8.5)
