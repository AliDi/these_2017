clear all
%close all

figure
semilogx(1,1)
hold on

clear all
load('../MCMC_sparse/data/MCMC_Mw_alphaVP.mat');
N=200; %nombre d'échantillon conservé après burnin
% 
% for i=1:length(Mw)
% 	for n = 1:N
% 		err_tmp(n)=norm( save_d_ref(:,i) - d_mcmc(:,end-N+n,i) ) /  norm( save_d_ref(:,i));
% 	end
% 	err_moy(i) = mean(err_tmp);
% 	stdd(i)=prctile(err_tmp,95)-err_moy(i);
% end

d_mcmc_moy = squeeze(mean(d_mcmc(:,end-200:end,:),2));
for i=1:length(Mw)
	stdd(:,i)=prctile(d_mcmc(:,end-N:end,i),95,2); %percentile suivant la deuxième dimension
    std_err(i) = norm( save_d_ref(:,i) - stdd(:,i) ) /  norm( save_d_ref(:,i));

end


%semilogx(Mw,10*log10(err_moy));
shadedErrorBar(Mw,10*log10(err_moy),10*log10(err_moy+stdd)-10*log10(err_moy))

%plot_fig(gcf,8,6)
%matlab2tikz('IC_Mw.tex')
