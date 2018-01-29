load('EM/EM_K.mat');
load('MCMC_sparse/MCMC_k.mat')
figure
hold on
plot(K,10*log10(err_EM(K)));



d_mcmc_moy = squeeze(mean(d_mcmc(:,end-200:end,:),2));
for i=1:length(K)
	err_moy(i) = norm( d_ref - d_mcmc_moy(:,i) ) /  norm( d_ref);
end
plot(K,10*log10(err_moy))

xlabel('Number of factors $\kappa$')
ylabel('Relative error on diag($\bo{S}_{aa}$) (dB)')
legend('EM','MCMC')
xlim([0 93])
ylim([-26 0])
plot_fig(gcf,8,6)
matlab2tikz('EM_MCMC_factors.tex')
