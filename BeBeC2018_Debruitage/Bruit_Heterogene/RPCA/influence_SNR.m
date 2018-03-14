%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =20;
Mw=10^4;
rho=0;
SNR = -10:10;
lambda=0:0.01:1;
for j=1:length(lambda)
    j
    for i=1:length(SNR)
        %%% Generate data
        [Sq Sy(:,:,i) Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw,SNR(i)-10);
        
        CSM = Sy(:,:,i);
        d_ref(:,i,j)=real(diag(Sp));
       
        %%%--------------------------------------------------------------------------------------------
        %%% RPCA solved with proximal gradient
        %%%--------------------------------------------------------------------------------------------
        [A(:,:,i,j) E(:,:,i,j) Nit out] = proximal_gradient_rpca(CSM , lambda(j), 300,1e-7,-1,-1,-1,-1);

        d_pg(:,i,j) = real( diag(A(:,:,i,j)) );	
        err_pg(i,j)= norm( (d_ref(:,i,j) - d_pg(:,i,j))) / norm(d_ref(:,i,j));
    end
end
save('RPCA_SNR','err_pg','lambda','d_pg','d_ref','SNR','A','E','Sy');


