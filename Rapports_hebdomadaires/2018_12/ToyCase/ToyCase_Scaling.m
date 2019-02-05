%Example from Rééchantillonnage de l'échelle dans les algorithmes MCMC pour les problèmes inverses bilinéaires, Veit 2008

clear all; %close all;
addpath('/home/adinsenmey/Bureau/these_2017/MCMC_Samplers/')
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')


Isnap = 10000; %nb d'échantillon
K=5; %nb de facteurs
M=20; % nb d'observations
SNR=-10 %en db

%{
graine=1;
load('s')
rng(s(graine));
L=(randn(M,K)+1i*randn(M,K))./sqrt(2*K);
gamma2=fliplr(1:K);
for i=1:length(gamma2)
    c(i,:)=(randn(Isnap,1) +1i*randn(Isnap,1)).*sqrt(gamma2(i)/2);
end
n=randn(M,Isnap) + 1i*randn(M,Isnap);
%}


load('L');
gamma2=[5 4 3 2 1];
load('c'); 
load('n');


%%
sig2 = (rms(L*c,2).*10.^((-SNR)./20)).^2;
n=n.*sqrt(sig2/2);

Sn=n*n'./Isnap;
Sc=c*c'./Isnap;
Sp=L*Sc*L';

y=L*c + n;
Sy=Sp+Sn;%y*y'./Isnap;
Sy=(Sy+Sy')/2;


%% Initialisation
K_est=M-1;
graine=1;
load('s')
rng(s(graine));
Ini.Lambda(:,:)=(randn(M,K_est)+1i*randn(M,K_est))./sqrt(2*K_est);
noise=real(mean(diag(Sy(:,:))));
%[a.sig2 , b.sig2] = Convert_InvGamma(noise,100*noise);
a.sig2 = 0.01;
b.sig2 = 0.01;
Ini.sig2=real(diag(Sy));
[a.beta2 , b.beta2] = Convert_InvGamma(1,10);
gamma2_mean =real(mean(diag(Sy(:,:))));
Ini.gamma2  = 0.3*gamma2_mean;
%[a.gamma2,b.gamma2] = Convert_InvGamma(0*gamma2_mean,100*gamma2_mean);   % hyperparamètres pour la variance des facteurs
a.gamma2 = 0.01;
b.gamma2=0.01;
[a.l b.l]=Convert_Beta(0.3,0.2);
a.l=1; b.l=1;

%% Inférence
Nrun=5000;
option.ref=0; option.noise='hetero'; option.marg='on'

% Pour la version sparse 3
for j=1:1
    %alpha2_mean(:,j) = svd(Sy(:,:,j+option.ref));%exp(-kappa*(0:K_est-1)'/(K_est-1));
    %alpha2_mean(:,j) = alpha2_mean(:,j)./max(alpha2_mean(:,j));
    a.alpha(:,j) = 2*ones(K_est,1);% 1./alpha2_mean(1:K_est,j);      % hyper-hyper-param�tres sur alpha
end
Ini.alpha(1,:,:) = 1./a.alpha;

%modele=2; [Sc_est,Lambda_est,q_est,l_est,beta2_est,gamma2_est, sig2_est,s] = MCMC_AnaFac_Quad_BernGauss_multiregime(Sy,K_est,a,b,Isnap,Nrun,option,Ini);
%modele=2; [Sc_est,Lambda_est,q_est,l_est,beta2_est,gamma2_est, sig2_est,s] = MCMC_AnaFac_Quad_BernGauss_multiregime_varL1surlq(Sy,K_est,a,b,Isnap,Nrun,option,Ini);
modele=1; [Sc_est,Lambda_est,q_est,beta2_est,gamma2_est, sig2_est] = MCMC_AnaFac_Quad_Sparse3_multiregime(Sy,K_est,a,b,Isnap,Nrun,option,Ini);
%%

for jj=1:Nrun
    tmp_Sp(:,:,jj)=squeeze(Lambda_est(jj,:,:))*diag(q_est(jj,:))*squeeze(Sc_est(jj,:,:))*diag(q_est(jj,:))*squeeze(Lambda_est(jj,:,:))';
    tmp_Sn(:,:,jj)=diag(beta2_est(jj).*sig2_est(jj,:));
end
Sp_est(:,:)=mean(tmp_Sp(:,:,round(Nrun/2):end),3);
Sn_est(:,:)=mean(tmp_Sn(:,:,round(Nrun/2):end),3);
Sy_est=Sp_est(:,:)+Sn_est(:,:);

err_Sp=norm(vec(Sp_est(:,:))-vec(Sp(:,:)))/norm(vec(Sp(:,:)));
err_reconstruction=norm(vec(Sy_est)-vec(Sy(:,:)))/norm(vec(Sy(:,:)));


%% Affichage
if modele==2
    figure
    subplot(231)
    plot(beta2_est.*mean(sig2_est,2))
    hold on
    plot(1:Nrun,ones(1,Nrun).*mean(diag(Sn)),'--k')
    xlabel("It\'erations")
    ylim(mean(diag(Sn))*[1-0.05 1+0.05])
    title('Valeur moyenne du bruit')

    subplot(232)
    plot(gamma2_est)
    title('Variance des facteurs : $\gamma^2$')

    subplot(234)
    plot(squeeze(sum(q_est,2)))
    hold on
    plot(1:Nrun,ones(1,Nrun).*K,'--k')
    title('Nombre de binaires non-nuls')

    subplot(235)
    plot(l_est)
    title('Param\`etre de parcimonie $l$','interpreter','latex')

    subplot(233)
    set(gca, 'visible', 'off')
    %title({['Erreur rms sur $S_p =$ ' num2str(10*log10(err_Sp)) ' dB'],''})
    text(0,0, {['Erreur rms sur'],[' $S_p =$ ' num2str(10*log10(err_Sp)) ' dB'],...
        ['graine al\''eatoire no ' num2str(graine)]})
    %set(findall(gca, 'type', 'text'), 'visible', 'on')

    subplot(236)
    for k=1:K_est
        for i=1:Nrun
            z(i,k)=norm(Lambda_est(i,:,k));
        end
    end
    plot(z)
    title('$||L_k||_2$')
    xlabel("It\'erations")
end



if modele==1
        figure
    subplot(231)
    plot(beta2_est.*mean(sig2_est,2))
    hold on
    plot(1:Nrun,ones(1,Nrun).*mean(diag(Sn)),'--k')
    xlabel("It\'erations")
    ylim(mean(diag(Sn))*[1-0.05 1+0.05])
    title('Valeur moyenne du bruit')

    subplot(232)
    plot(gamma2_est)
    title('Variance des facteurs : $\gamma^2$')

    subplot(234)
    plot(squeeze(q_est))
    title('q')
    
    subplot(233)
    set(gca, 'visible', 'off')
    %title({['Erreur rms sur $S_p =$ ' num2str(10*log10(err_Sp)) ' dB'],''})
    text(0,0, {['Erreur rms sur'],[' $S_p =$ ' num2str(10*log10(err_Sp)) ' dB'],...
        ['graine al\''eatoire no ' num2str(graine)]})
    set(findall(gca, 'type', 'text'), 'visible', 'on')

    subplot(236)
    for k=1:K_est
        for i=1:Nrun
            z(i,k)=norm(Lambda_est(i,:,k));
        end
    end
    plot(z)
    title('$||L_k||_2$')
    xlabel("It\'erations")
    
    subplot(235)
    set(gca, 'visible', 'off')
end

for i=1:6
    subplot(2,3,i)
    plot_fig(gcf,20,12)
end


%{
figure
plot(beta2_est.*mean(sig2_est,2))
hold on
plot(Nrun,mean(diag(Sn)),'*')
plot(1:Nrun,ones(1,Nrun).*mean(diag(Sn)))

figure
plot(squeeze(sum(q_est,2)))

figure
subplot(1,2,2)
imagesc(squeeze(real(Sp_est)))
subplot(1,2,1)
imagesc(squeeze(real(Sp)))

%}








