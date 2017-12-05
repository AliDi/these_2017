%comparison of several denoising algorithms in term of  : 
%-error : e = | diag_reconstructed - diag(Sp) | / |diag(Sp)|
%-convergence : c = | diag_rec(k-1) - diag(Sp)| / |diag_rec(k) - diag(Sp)|
%%
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 3800;
Mw=9000;
Nsrc = 5;
rho=0;
SNR=-20:1:20;

for i=1:length(SNR)
i
    %%% Generate data
    [Sq b Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw);
    Sy(:,:,i)=Sp+diag(diag((Sn)));

    CSM = Sy(:,:,i);
    d_ref(:,i)=real(diag(Sp));

    %%%--------------------------------------------------------------------------------------------
    %%% SLDR solved with proximal gradient
    %%%--------------------------------------------------------------------------------------------
	lambda= 0.5;
	[A E Nit out A_all] = proximal_gradient_rpca(CSM , lambda, 300,1e-8,-1,-1,-1,-1);
	
	d_pg(:,i) = real( diag(A) );	
	err_pg(i)= norm( diag(d_ref(:,i) - d_pg(:,i))) / norm(d_ref(:,i));

end

save('RPCA_SNR_5src','err_pg','SNR','d_pg','d_ref');

%figure
%plot(SNR,10*log10(err_pg))
%xlabel('SNR (dB)')
%ylabel('Relative error on diag($\boldmath{S_y}$) (dB)')
%xlim([-20 20])
%plot_fig(gcf,11,8.5)

