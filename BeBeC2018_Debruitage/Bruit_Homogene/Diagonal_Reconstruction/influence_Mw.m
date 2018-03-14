%Alice Dinsenmeyer; hiver 2017-2018
clear all
%close all

addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),50));
rho=0;
SNR = 10;
for j=1:length(Mw)
    j
	for i=1:length(Nsrc)
	i
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j),SNR);
		CSM = Sy(:,:);
		d_ref(:,i,j)=real(diag(Sp));

		%%%--------------------------------------------------------------------------------------------
		%%% Hald
		%%%--------------------------------------------------------------------------------------------
		cvx_quiet('true')
		cvx_precision('high')
		[CSM_cvx d1 cvx_it]=CSMRecHald(CSM);     
		d_cvx(:,i,j)=diag(CSM_cvx);
		err_cvx(i,j)= norm( d_ref(:,i,j)-real(d_cvx(:,i,j)) ,2)/norm(d_ref(:,i,j),2); 

		%%%--------------------------------------------------------------------------------------------
		%%% Alternating projections
		%%%--------------------------------------------------------------------------------------------
	    [CSM_it, n, errec, D_AP] = recdiag(CSM,1,1000,1e-7,30); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
	    d_it(:,i,j) = real(diag(CSM_it));
	    err_it(i,j) = norm( d_ref(:,i,j)-real(d_it(:,i,j))) / norm(d_ref(:,i,j));
	end
end
save('cvx_Mw','err_cvx','Mw','Nsrc','d_cvx', 'd_ref');
save('AP_it_Mw','err_it','Mw','Nsrc','d_it','d_ref');

%%%--------------------------------------------------------------------------------------------
%%% Dougherty
%%%--------------------------------------------------------------------------------------------
% Moins de points de calcul car algo plus long
Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),25));

for j=1:length(Mw)
    j
	for i=1:length(Nsrc)
	i
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j));
		CSM = Sy(:,:);
		d_ref(:,i,j)=real(diag(Sp));
		
		 %linprog (Dougherty)
	     [CSM_linprog ii] = recdiagd(double(CSM),500,60);
	     d_linprog(:,i,j) = real(diag(CSM_linprog));
	     err_linprog(i,j)= norm( d_ref(:,i,j)-real(d_linprog(:,i,j))) / norm(d_ref(:,i,j)) ;
	end
end
save('linprog_Mw','err_linprog','Mw','Nsrc','d_linprog','d_ref');
