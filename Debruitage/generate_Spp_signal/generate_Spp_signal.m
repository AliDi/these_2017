% Alice Dinsenmeyer, 24/10/17
%Generates cross spectral matrix depending on : 
%-source positions, amplitude
%-antenna location
%-propagator G (function of frequency)
%-number of samples and snapshots%

clear all; 
close all;

%rng(s)
%%%--------------------------------------------------------------------------------------------
%%% Acquisition parameters
%%%--------------------------------------------------------------------------------------------

%time window parameters
Mw=2^11; %number of time windows 
overlap = 0; % overlaping of time windows : 0<overlap<1
Mw=Mw-overlap*(Mw-1); %effective number of time windows

%Frequencies of interest
freqs=1000;%[ 100 1000 3000];
Nfreq = length(freqs);

fftsign=+1;

%%%--------------------------------------------------------------------------------------------
%%% Antenna parameters
%%%--------------------------------------------------------------------------------------------

Nmic = 196;

z_mic= 0.8 * ones(1,Nmic);
x_mic = linspace(0,0.6,sqrt(Nmic));
y_mic = linspace(-.3,.3,sqrt(Nmic));

[x_mic y_mic] = meshgrid(x_mic,y_mic);

r_mic = [x_mic(:) y_mic(:) z_mic(:)];

%%%--------------------------------------------------------------------------------------------
%%% Source parameters
%%%--------------------------------------------------------------------------------------------

Nsrc = 1; %number of sources

%positions (used for propagation)
z_src = 1 * ones(1,Nsrc); 
x_src =logspace(log10(0.10),log10(0.60),Nsrc); 
y_src = zeros(1,Nsrc);

r_src = [x_src' y_src' z_src']; %coordinates of real sources

%strength coefficients
qrms = 1*ones(Nsrc,1); %rms value supposed to be constant on frequency

%%%--------------------------------------------------------------------------------------------
%%% Noise parameters
%%%--------------------------------------------------------------------------------------------

SNR = -20; %signal to noise ration in dB

%%%--------------------------------------------------------------------------------------------
%%% Calculate noisy spectra 
%%%--------------------------------------------------------------------------------------------

Q(:,:,:)=qrms.*(randn(Nsrc,Mw,Nfreq) + 1i*randn(Nsrc,Mw,Nfreq))/sqrt(2); %TF of q(t) for Mw snapshots and Nfreq frequencies

for f=1:Nfreq

	%choose a propagator
	%G_src_mic(:,:,f) = randn(Nmic,Nsrc)+1i*randn(Nmic,Nsrc); %random mixing matrix
    %G_src_mic(:,:,f) = GreenFreeFieldUniformFlow(r_src , r_mic , freqs(f) , -fftsign ,Mach); %monopole in uniform flow
    G_src_mic(:,:,f) = GreenFreeField(r_src, r_mic , freqs(f) , -fftsign); %monopole in free field
        
    %Calculate propagated signals
	P(:,:,f) = G_src_mic(:,:,f)*Q(:,:,f);
	
	%Generate Noise
	nrms = mean(abs(P(:,:,f))).*10.^((-SNR)./20); %rms value
	
	N(:,:,f) = nrms.*(randn(Nmic,Mw) + 1i*randn(Nmic,Mw))/sqrt(2);
	
	%Add noise to signals 	
	Y(:,:,f) = P(:,:,f) + N(:,:,f);
	
	%Calculate mean spectra 
	Sy(:,:,f) = Y(:,:,f) * Y(:,:,f)' /Mw; %signal + noise
	Sp(:,:,f) = P(:,:,f) * P(:,:,f)' /Mw; %mic signals
	Sn(:,:,f) = N(:,:,f) * N(:,:,f)' /Mw; %noise only	 
end

%%%--------------------------------------------------------------------------------------------
%%% Display data
%%%--------------------------------------------------------------------------------------------

figure
plot3(r_mic(:,1),r_mic(:,2),r_mic(:,3),'ob')
hold on
plot3(x_src,y_src,z_src,'or')
title('Configuration')
legend('mic positions','src positions')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

%%%--------------------------------------------------------------------------------------------
%%% Export data
%%%--------------------------------------------------------------------------------------------
Sp_filename='Sp.h5';

if exist(Sp_filename)
	d=input(['Delete ' Sp_filename ' (y/n) ? '],'s');
	
	if d=='y'
		delete(Sp_filename);
		h5create(Sp_filename,'/Spp',size(Sp));
		h5write(Sp_filename,'/Spp',Sp);

		h5create(Sp_filename,'/frequencies',size(freqs));
		h5write(Sp_filename,'/frequencies',freqs);
	else 
		disp('Data not saved');
	end
end









%A supprimer
%====Uncorrelated Gaussian signals
%Sqq=zeros(Nsrc,Nsrc,Nfreq);
%for f=1:Nfreq
%	%cross-spectra
%	for i = 1:Nsrc
%		for j = 1:i-1
%		    Sqq(i,j,f)= qrms(i)*qrms(j)*(randn+ 1i * randn) / ( sqrt(Mw)*sqrt(2) ) ;
%		    Sqq(j,i,f)= conj(Sqq(i,j,f));
%		end
%	end
%	
%	%add auto-spectra
%	Sqq(:,:,f) = Sqq(:,:,f) + diag(qrms.^2 .* ( 1 + randn(Nsrc,1)/sqrt(Mw)) );
%end



%%%--------------------------------------------------------------------------------------------
%%% Generate noise matrix Sn
%%%--------------------------------------------------------------------------------------------
%SNR is the same for all frequencies and vary around +- 1dB over mics 


%for f=1:Nfreq
%	%noise level on each mic
%	 %rms value
%	 nrms = mean(real(diag(Spp(:,:,f)))).*10.^((-SNR+randn(Nmic,1))./10); %rms value
%	
%	%cross-spectra
%	for i = 1:Nmic
%		for j = 1:i-1
%		    Snn(i,j,f)= nrms(i)*nrms(j)*(randn+ 1i * randn) / ( sqrt(Mw)*sqrt(2) ) ;
%		    Snn(j,i,f)= conj(Snn(i,j,f));
%		end
%	end
%	
%	%add auto-spectra
%	Snn(:,:,f) = Snn(:,:,f) + diag(nrms.^2 .* ( 1 + randn(Nmic,1)/sqrt(Mw)) );
%end
%====    %Spp(:,:,f) = G_src_mic(:,:,f)*Sqq(:,:,f)*G_src_mic(:,:,f)';
%%%--------------------------------------------------------------------------------------------
%%% Add noise to signal 
%%%--------------------------------------------------------------------------------------------
%Sy = Snn + Spp;

