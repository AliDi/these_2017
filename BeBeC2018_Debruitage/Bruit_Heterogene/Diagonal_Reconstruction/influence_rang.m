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
Mw=10^4;
rho=0;
SNR = 10;

for j=1:length(Mw)
    j
	for i=1:length(Nsrc)
	i
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j),SNR-10);
		CSM = Sy(:,:);
		%d_ref(:,i,j)=real(diag(Sy-Sn));
		d_ref(:,i,j)=real(diag(Sp));
		%%%--------------------------------------------------------------------------------------------
		%%% SLDR solved with proximal gradient
		%%%--------------------------------------------------------------------------------------------
		%lambda= 0.1;
		%[A E Nit out A_all] = proximal_gradient_rpca(CSM , lambda, 70,1e-10,-1,-1,-1,-1);

		%d_pg = real( diag(A_all(:,:,end)));	
		%err_pg(i)= norm( diag(d_ref - d_pg) ,1) / norm(d_ref,1);

		%%%--------------------------------------------------------------------------------------------
		%%% Alternating projections
		%%%--------------------------------------------------------------------------------------------
		%cvx (Hald)
		cvx_quiet('true')
		cvx_precision('high')
		[CSM_cvx d1 cvx_it]=CSMRecHald(CSM);     
		d_cvx(:,i,j)=diag(CSM_cvx);
		err_cvx(i,j)= norm( d_ref(:,i,j)-real(d_cvx(:,i,j)) ,2)/norm(d_ref(:,i,j),2); 


 		 %linprog (Dougherty)
% 	     [CSM_linprog ii] = recdiagd(double(CSM),500,60);
% 	     d_linprog(:,i,j) = real(diag(CSM_linprog));
% 	     err_linprog(i,j)= norm( d_ref(:,i,j)-real(d_linprog(:,i,j))) / norm(d_ref(:,i,j)) ;

		%Alternating projection
 	    %[CSM_it, n, errec, D_AP] = recdiag(CSM,1,1000,1e-7,30); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
 	    %d_it(:,i,j) = real(diag(CSM_it));
 	    %err_it(i,j) = norm( d_ref(:,i,j)-real(d_it(:,i,j))) / norm(d_ref(:,i,j));

		%Version Jérôme (remove average value of the K greatest SV)
		%[Sn,L] =
	%         SS_CSM_Fit(CSM,92);
		%d_jerome(:,i) = real(diag(L*L'));
		%d_jerome2(:,i) = real(diag(CSM-diag(Sn)));
		%err_jerome(i) = norm( d_ref(:,i)-real(d_jerome(:,i))) / norm(d_ref(:,i));

	end
end
save('cvx_rang','err_cvx','Mw','Nsrc','d_cvx', 'd_ref');
% 
% save('AP_it_rang','err_it','Mw','Nsrc','d_it','d_ref');
% 
% save('linprog_rang','err_linprog','Rang2','d_linprog','d_ref');
% 


