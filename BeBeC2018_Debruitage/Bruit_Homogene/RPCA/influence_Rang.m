%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =1:93;
Mw=10841;%10^4;
rho=0;
SNR = 10;
lambda=0:0.01:1;

for i=1:length(Nsrc)
i
    %%% Generate data
    [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw , SNR );
    
    CSM = Sy(:,:);
    d_ref(:,i)=real(diag(Sp));
    %%%--------------------------------------------------------------------------------------------
    %%% RPCA solved with proximal gradient
    %%%--------------------------------------------------------------------------------------------
    for j=1:length(lambda)
    
    	[A E Nit out] = proximal_gradient_rpca(CSM , lambda(j), 300,1e-7,-1,-1,-1,-1);

    	d_pg(:,i,j) = real( diag(A) );	
    	err_pg(i,j)= norm( (d_ref(:,i) - d_pg(:,i,j))) / norm(d_ref(:,i));
    end
end

for i=1:length(Nsrc)
	[err_pg_opt(i) b(i)]=min(err_pg(i,:));
	
end

save('RPCA_rang','err_pg','lambda','d_pg','d_ref','err_pg_opt');


