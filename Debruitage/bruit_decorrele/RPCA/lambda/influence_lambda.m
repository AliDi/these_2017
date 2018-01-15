%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 3800;
Nsrc =50;
Mw=9000;
rho=0;
SNR = 10;

lambda=0.001:0.01:1;

%%% Generate data
[Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw);
Sy(:,:)=Sp+diag(diag((Sn)));
CSM = Sy(:,:);
d_ref=real(diag(Sp));

for i=1:length(lambda)
i

	%%%--------------------------------------------------------------------------------------------
	%%% RPCA solved with proximal gradient
	%%%--------------------------------------------------------------------------------------------
	%lambda= 0.5;
	[A E Nit out] = proximal_gradient_rpca(CSM , lambda(i), 70,1e-10,-1,-1,-1,-1);

	d_pg(:,i) = real( diag(A) );	
	err_pg(i)= norm( diag(d_ref - d_pg(:,i))) / norm(d_ref);

end

save('RPCA_lambda','err_pg','lambda','d_pg','d_ref');



