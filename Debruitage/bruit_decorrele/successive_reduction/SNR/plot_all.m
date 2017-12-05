cvx93=load('cvx_Rang93.mat');
cvx5=load('cvx_Rang5.mat');

ap5=load('AP_it_Rang5.mat');
ap93=load('AP_it_Rang93.mat');

lp5=load('linprog_snr_5src.mat');
lp93=load('linprog_snr_150src.mat');


figure
hold on
plot(cvx93.SNR,10*log10(cvx93.err_cvx))
plot(ap93.SNR,10*log10(ap93.err_it))
plot(lp93.SNR,10*log10(lp93.err_linprog))


plot(cvx5.SNR,10*log10(cvx5.err_cvx))
plot(ap5.SNR,10*log10(ap5.err_it))
plot(lp5.SNR,10*log10(lp5.err_linprog))

ylabel('Relative error on diag($\boldmath{S_y}$) (dB)')
xlabel('SNR (dB)')

leg=legend('Hald, $R=93$','Basic AP, $R=93$','Dougherty, $R=93$','Hald, $R=5$','Basic AP, $R=5$','Dougherty, $R=5$','location','southeastoutside');


plot_fig(gcf,500,300)