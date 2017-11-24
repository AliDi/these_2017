function [l cdf]=distrib_vp_noise(q,p,s)

%Estimation de la répartition des valeurs propres d'une matrice de covariance de bruit gaussien
%/!\ le nombre de snapshots doit être strictement supérieur au nombre de microphones
%voir DISTRIBUTION OF EIGENVALUES FOR SOME SETS OF RANDOM MATRICES, Marchenko & Pastur 1967
%et Eigenvalues of the sample covariance matrix for a towed array, Gerstoft et al.2012

%q : nombre de microhpones
%p : nombre de snapshots
%s : noise variance

nu=q/p;

%%% Calcul des PDF et CDF (DISTRIBUTION OF EIGENVALUES FOR SOME SETS OF RANDOM MATRICES, Marchenko & Pastur 1967)
n=10000;  %nombre de point pour le calcul de la pdf et de la cdf
lm=s*(1-sqrt(nu))^2; %borne inférieure des VP
lp=s*(1+sqrt(nu))^2; %borne supérieure des VP
l=linspace(lm,lp,n); % VP

%pdf
P=sqrt((lp-l).*(l-lm))./(2*pi*nu.*l*s); %pdf

%cdf
for i=1:n
	cdf(i)=sum(P(1:i));
end
cdf=cdf/max(cdf); %normalisation
cdf=cdf*(q-1)+1;

end
