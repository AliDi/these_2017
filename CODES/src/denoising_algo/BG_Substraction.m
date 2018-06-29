function [Gsig]=BG_Substraction(G, Gd)

%Background substraction from
%Advanced Background Subtraction Applied to Aeroacoustic Wind Tunnel Testing, de C. J. Bahr et W. C. Horne
%https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005974.pdf
%
% G = Gsig + Gd
%
%INPUTS : 
%G : (MxM) cross-spectral matrix (CSM) to be denoised
%Gd : (MxM) background noise CSM


%% Check if Gd is Hermitian positive semidefinite
[~ , test_pos]=chol(Gd);
if test_pos %if Gd negative
	Gd=nearestSPD(Gd);
end

%% Eigendecomposition of Gd
[Xd, Ld]=eig(Gd);

%% Prewhitening operator Bd
Bd = Xd*diag(diag(Ld).^(-1/2));

%% Prewhiten G
Ghat=Bd'*G*Bd;

%% Eigendecomposition of Ghat
[Xhat, Lhat]=eig(Ghat);

%% Extract positive eigenvalues for signal
Lsighat = diag(Lhat-1);
Lsighat(Lsighat<0) = 0;
Lsighat=diag(Lsighat);


%% Compute Gsighat
Gsighat=Xhat * Lsighat * Xhat';

%% Inverse prewhitening operator
Bdinv=inv(Bd);
Gsig = Bdinv' * Gsighat * Bdinv;

