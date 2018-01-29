clear all
close all

figure

% alternating projections
load('cvx_Mw.mat')
semilogx(Mw,10*log10(err_cvx(1,:)))

hold on


%RPCA
load('RPCA_Mw.mat')
semilogx(Mw,10*log10(err_pg_opt))
semilogx(Mw,10*log10(err_pg(:,11)))

%MCMC
load('MCMC_Mw.mat');
d_mcmc_moy = squeeze(mean(d_mcmc(:,end-200:end,:),2));
for i=1:length(Mw)
	err_moy(i) = norm( d_ref(:,i) - d_mcmc_moy(:,i) ) /  norm( d_ref(:,i));
end
semilogx(Mw,10*log10(err_moy));


%EM
%load('EM/EM_Mw.mat');
%semilogx(Mw,10*log10(err_EM(1,:)));

xlabel('Number of snapshots $N_s$')
ylabel('Relative error on diag($\bm{S}_{aa}$) (dB)')
legend('Hald','RPCA, $\lambda_{opt}$','RPCA, $\lambda=1/\sqrt{M}$','MCMC')

%plot_fig(gcf,8,6)
%matlab2tikz('all_Mw.tex')
