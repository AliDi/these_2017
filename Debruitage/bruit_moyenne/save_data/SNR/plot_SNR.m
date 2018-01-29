clear all
close all

figure
hold on

% alternating projections
load('cvx_SNR.mat')
plot(SNR,10*log10(err_cvx))

%load('EM/EM_Mw.mat');
%plot(10*log10(err_EM))

%RPCA
load('RPCA_SNR.mat')

for i=1:length(SNR)
	err_pg_opt(i)=min(err_pg(i,:));	
end
plot(SNR,10*log10(err_pg_opt))
plot(SNR,10*log10(err_pg(:,11))) %Wright's parameter'

%MCMC
load('MCMC_SNR.mat');
d_mcmc_moy = squeeze(mean(d_mcmc(:,end-200:end,:),2));
for i=1:length(SNR)
	err_moy(i) = norm( d_ref(:,i) - d_mcmc_moy(:,i) ) /  norm( d_ref(:,i));
end
plot(SNR,10*log10(err_moy));

xlabel('SNR (dB)')
ylabel('Relative error on diag($\bo{S}_{aa}$) (dB)')
legend('Hald','RPCA, $\lambda_{opt}$','RPCA, $\lambda=1/\sqrt{M}$','MCMC')

%plot_fig(gcf,8,6)
%matlab2tikz('all_SNR.tex')
