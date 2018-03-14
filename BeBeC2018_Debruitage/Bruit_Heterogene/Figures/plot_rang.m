clear all
%close all

figure
hold on

% alternating projections
clear all
load('../Diagonal_Reconstruction/data/Rang/cvx_rang.mat')
plot(10*log10(err_cvx))

%EM
%load('EM/EM_Mw.mat');
%plot(10*log10(err_EM))

%RPCA
clear all
load('../RPCA/data/RPCA_rang.mat')

for i=1:93
	err_pg_opt(i)=min(err_pg(i,:));	
end
plot(10*log10(err_pg_opt))
plot(10*log10(err_pg(:,11)))

%MCMC
clear all
 load('../MCMC_sparse/data/MCMC_Nsrc_alphaVP.mat');
 d_mcmc_moy = squeeze(mean(d_mcmc(:,end-200:end,:),2));
 for i=1:93
 	err_moy(i) = norm( save_d_ref(:,i) - d_mcmc_moy(:,i) ) /  norm( save_d_ref(:,i));
 end
 plot(10*log10(err_moy));

xlabel('Rank of $\bm{S}_{aa}$')
ylabel('Relative error on diag($\bm{S}_{aa}$) (dB)')
legend('Hald','RPCA, $\lambda_{opt}$','RPCA, $\lambda=1/\sqrt{M}$','MCMC')

%plot_fig(gcf,8,6)
%matlab2tikz('all_rang.tex')
