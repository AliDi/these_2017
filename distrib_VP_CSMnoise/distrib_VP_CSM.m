%Estimation de la répartition des valeurs propres d'une matrice de covariance de bruit gaussien
%/!\ le nombre de snapshots doit être strictement supérieur au nombre de microphones
%voir DISTRIBUTION OF EIGENVALUES FOR SOME SETS OF RANDOM MATRICES, Marchenko & Pastur 1967
%et Eigenvalues of the sample covariance matrix for a towed array, Gerstoft et al.2012


clear all;
%close all;
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
format long	

%%% Paramètres à faire varier
q=196; %nombre de microphones
p= 196 * 10; %nombre de snapshots
nu=q/p;

%%% Calcul des PDF et CDF (DISTRIBUTION OF EIGENVALUES FOR SOME SETS OF RANDOM MATRICES, Marchenko & Pastur 1967)
n=10000;  %nombre de point pour le calcul de la pdf et de la cdf
lm=(1-sqrt(nu))^2; %borne inférieure des VP
lp=(1+sqrt(nu))^2; %borne supérieure des VP
l=linspace(lm,lp,n); % VP

%pdf
P=sqrt((lp-l).*(l-lm))./(2*pi*nu.*l); %pdf

%cdf
for i=1:n
	cdf(i)=sum(P(1:i));
end
cdf=cdf/max(cdf); %normalisation

figure(1);
plot(l,P);
title('Densit\''e de probabilit\''e des valeurs propres')
ylabel('$p(\lambda)$')
xlabel('$\lambda$')

%%% Comparaison avec un résultat expériemental
N(:,:) = (randn(q,p) + 1i*randn(q,p))/sqrt(2); %bruit 
Sn(:,:) = N(:,:) * N(:,:)' /(p); %matrice de covariance de bruit

figure(2);
hold on 
plot(cdf*(q-1)+1,10*log10(l)) %CDF normalisée
plot(sort(10*log10(real(eig(Sn)))));
ylabel('Amplitude des valeurs propres en dB')
legend('Pr\''ediction','Exp\''erience')
title(['Nombre de snapshots : ' num2str(1/nu) '$\times N_{mic}$'])




