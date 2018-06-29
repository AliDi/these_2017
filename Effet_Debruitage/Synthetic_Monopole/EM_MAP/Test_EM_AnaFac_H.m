% Teste le code EM_AnaFac_H.m
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/AnaFac_H/')
clear all
close all
%load('s')
%rng(s)

% Param�tres de simulation
M = 20;   % nombre de capteurs
N = 21;          % nombre de fonctions de base
K = 5;          % nombre de "sources" (nombre de facteurs)
Isnap = 1000;    % nombre de snapshots

% Fonctions de transfert
H = (randn(M,N) + 1i*randn(M,N))/sqrt(2);
%H = H/K;
%H = eye(M,N);

% Fonctions de base
Lambda = (randn(N,K) + 1i*randn(N,K))/sqrt(2);

% Mod�le direct
alpha2 = 1;
beta2 = 3;%1e-1;

c = sqrt(alpha2)*(randn(K,Isnap) + 1i*randn(K,Isnap))/sqrt(2);
n = sqrt(beta2)*(randn(M,Isnap) + 1i*randn(M,Isnap))/sqrt(2);
q = Lambda*c;
y = H*q + n;

% Matrice spectrale observ�e
Sy = y*y'/Isnap;

% Matrice spectrale sources
Sq = q*q'/Isnap;

% Matrice spectrale du bruit
Sn = n*n'/Isnap;

%Pour un bruit parfaitement décorrélé
%Sy=H*q*(H*q)'/Isnap+diag(diag(Sn));

figure
subplot(221),imagesc(real(Sq)),colorbar,title('Sq')
subplot(222),imagesc(real(H*Sq*H')),colorbar,title('H*Sq*H''')
subplot(223),imagesc(real(n*n'/Isnap)),colorbar,title('Sn')
subplot(224),imagesc(real(Sy)),colorbar,title('Sy')

%% EM 
%%%%%%%
option.max=1000; %max number of iterations
option.rerr=1e-3;

Ini.Syc= mean(diag(Sy))*ones(M,1);%beta2*ones(M,1); %bruit moyen
Ini.L=randn(N,K)+1i*randn(N,K); %fonctions de base

[L,sig2,beta2,flag,Sx,loglik] = EM_AnaFac_H(Sy,Isnap,K,H,option,Ini);

%%

figure
subplot(221),imagesc(real(L*Sx*L')),colorbar,title('L*Sx*L''')
subplot(222),imagesc(real(H*L*Sx*L'*H')),colorbar,title('H*L*Sx*L''*H''')
subplot(223),imagesc(real(diag(sig2))),colorbar,title('diag(sig2)')
subplot(224),imagesc(real(H*L*Sx*L'*H')+diag(sig2)),colorbar,title('H*L*Sx*L''*H'')+diag(sig2)')

figure;
subplot(311); hold on
plot(real(diag(Sq)),'-o'), plot(real(diag(L*Sx*L')))
title('Debits des sources')
%legend('Reels','Reconstruits')

subplot(312); hold on
plot(real(diag(Sn)),'o-'),plot(real((sig2)))
title('Bruit')
%legend('Reel','Reconstruit')

subplot(313); hold on
plot(real(diag(Sy)),'o-'),plot(real(diag(H*L*Sx*L'*H'))+(sig2))
title('Pression mesuree')
legend('Reelle','Reconstruite','Location','bestoutside')


