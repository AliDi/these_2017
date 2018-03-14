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
Mw=round(logspace(log10(10),log10(50000),40));
rho=0;
SNR = 10;

j=1;
for i=1:length(Mw)
i
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw(i) , SNR );
	
	CSM = Sy(:,:);
	d_ref(:,i,j)=real(diag(Sp));
			
	if Mw(i)<93
   		 k=Mw(i)-1;
	else
    	k=92;
	end
	k
	
	[Sn,L,D] = SS_CSM_Fit(CSM,k);

	d_SI(:,i)=real(diag(L*L'));
	err_SI(i) = norm( d_ref(:,i) - d_SI(:,i) ) /  norm( d_ref(:,i));



end
save('test_SI_Mw','err_SI','Mw','Nsrc','d_SI', 'd_ref');
figure
semilogx(Mw,10*log10(err_SI))



