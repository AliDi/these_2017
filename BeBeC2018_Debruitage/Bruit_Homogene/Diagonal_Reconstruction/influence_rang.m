%Alice Dinsenmeyer; hiver 2017-2018
%%
clear all
%close all

addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 15000;
Nsrc =1:1:93;
Mw=10^4;
rho=0;
SNR = 10;
j=1; %variable inutile
for i=1:length(Nsrc)
	i
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j),SNR);
	Rang2(i,j)=rank(Sp);
	CSM = Sy(:,:);
	d_ref(:,i,j)=real(diag(Sp));
	
	%%%--------------------------------------------------------------------------------------------
	%%% Hald
	%%%--------------------------------------------------------------------------------------------
	cvx_quiet('true')
	cvx_precision('high')
	[CSM_cvx d1 cvx_it]=CSMRecHald(CSM);     
	d_cvx(:,i,j)=real(diag(CSM_cvx));
	err_cvx(i,j)= norm( d_ref(:,i,j)-real(d_cvx(:,i,j)) ,2)/norm(d_ref(:,i,j),2); 

	%%%--------------------------------------------------------------------------------------------
	%%% Alternating projections
	%%%--------------------------------------------------------------------------------------------
    [CSM_it, n, errec, D_AP] = recdiag(CSM,1,1000,1e-7,30); %enlever les VP negatives puis recalculer les VP en interchangeant la diagonale
    d_it(:,i,j) = real(diag(CSM_it));
    err_it(i,j) = norm( d_ref(:,i,j)-real(d_it(:,i,j))) / norm(d_ref(:,i,j));
end

save('cvx_rang','err_cvx','Mw','Nsrc','d_cvx', 'd_ref');
save('AP_it_rang','err_it','Mw','Nsrc','d_it','d_ref');

%%%--------------------------------------------------------------------------------------------
%%% Dougherty
%%%--------------------------------------------------------------------------------------------
Nsrc =1:2:96;
for i=1:length(Nsrc)
	i
	%%% Generate data
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j));
	Rang2(i,j)=rank(Sp);
	CSM = Sy(:,:);
	d_ref(:,i,j)=real(diag(Sp));
	
	 %linprog (Dougherty)
     [CSM_linprog ii] = recdiagd(double(CSM),500,60);
     d_linprog(:,i,j) = real(diag(CSM_linprog));
     err_linprog(i,j)= norm( d_ref(:,i,j)-real(d_linprog(:,i,j))) / norm(d_ref(:,i,j)) ;
end
save('linprog_rang','err_linprog','Rang2','d_linprog','d_ref');



