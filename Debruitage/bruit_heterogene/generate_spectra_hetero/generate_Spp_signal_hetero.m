function [Sq Sy Sp Sn] = generate_Spp_signal_hetero(freqs , Nsrc , rho , SNR, Mw, extra_SNR)
%function [Sq Sy Sp Sn] = generate_Spp_signal(freqs , Nsrc , rho , SNR)
%freqs : vector, frequencies
%Nsrc : number of sources
%rho : amount of source correlation. 0<rho<1 : partially correlated ; 0 : uncorrelated , ...
%SNR : SNR in dB added to all the microphones
%extra_SNR : SNR in dB added to 10 random mic
%
%
%outputs : 
%Sq : sources spectra
%Sy : signal + noise spectra
%Sp : signal
%Sn : noise

%Generates cross spectral matrix depending on : 
%-source positions, amplitude
%-antenna location
%-propagator G (function of frequency)
%-number of samples and snapshots%


%load('s')
%rng(s)
%rng(s)
%%%--------------------------------------------------------------------------------------------
%%% Acquisition parameters
%%%--------------------------------------------------------------------------------------------

%%time window parameters
%Mw=100*93; %number of time windows 
%overlap = 0; % overlaping of time windows : 0<overlap<1
%Mw=Mw-overlap*(Mw-1); %effective number of time windows

%Frequencies of interest
Nfreq = length(freqs);

fftsign=+1;

%%%--------------------------------------------------------------------------------------------
%%% Antenna parameters
%%%--------------------------------------------------------------------------------------------
%[x_mic y_mic] = meshgrid(x_mic,y_mic);
%r_mic = [x_mic(:) y_mic(:) z_mic(:)];

r_mic=load('spiral_array_coord.mat');
r_mic=r_mic.mic_info;

Nmic = size(r_mic,1);

%%%--------------------------------------------------------------------------------------------
%%% Source parameters
%%%--------------------------------------------------------------------------------------------
%positions (used for propagation)
y_src = 0;

%quadrillage
%z_src = linspace(-1.2,1.2,ceil(sqrt(Nsrc)));
%x_src =linspace(-0.5,0.5,floor(sqrt(Nsrc)));
%[x_src y_src z_src]=meshgrid(x_src,y_src,z_src);

%three lines
%z_src = [linspace(-1.2,1.2,floor(Nsrc/3)) linspace(-1.2,1.2,round(Nsrc/3))  linspace(-1.2,1.2,ceil(Nsrc/3))];
%x_src =[-0.5*ones(1,floor(Nsrc/3)) 0*ones(1,round(Nsrc/3)) 0.5*ones(1,ceil(Nsrc/3))]; 
%[x_src y_src]=meshgrid(x_src,y_src);

%two lines
%z_src = [linspace(-1.2,1.2,floor(Nsrc/2))  linspace(-1.2,1.2,ceil(Nsrc/2))];
%size(z_src)
%x_src =[-0.5*ones(1,floor(Nsrc/2)) 0.5*ones(1,ceil(Nsrc/2))]; 
%[x_src y_src]=meshgrid(x_src,y_src);

%one line : 
%z_src = linspace(-1.2,1.2,Nsrc);
%x_src =0.2;
%[x_src y_src z_src]=meshgrid(x_src,y_src,z_src);

%one line with rotation theta
theta=0.0175;
R=1.2;
z_src = linspace(-R,R,Nsrc);
x_src =0.2 + z_src*sin(theta);
z_src=z_src-z_src*(1-cos(theta));

[x_src y_src]=meshgrid(x_src,y_src);

r_src = [x_src(:) y_src(:) z_src(:)]; %coordinates of real sources

%strength coefficients
qrms = 1*ones(Nsrc,1); %rms value supposed to be constant on frequency

%source correlation (do not depend on frequency)
%rho : amount of source correlation. 0<rho<1 : partially correlated ; 0 : uncorrelated , ...
Qcorr=eye(Nsrc);
Qcorr(Qcorr==0)=rho;
Qcorr=chol(Qcorr)';

%%%--------------------------------------------------------------------------------------------
%%% Select randomly 10 noisy mic
%%%--------------------------------------------------------------------------------------------
nb_noisy_mic=10;
tmp=randperm(Nmic);
noisy_mic=tmp(1:nb_noisy_mic);



%%%--------------------------------------------------------------------------------------------
%%% Calculate noisy spectra 
%%%--------------------------------------------------------------------------------------------
Q(:,:,:)=qrms.*(randn(Nsrc,Mw,Nfreq) + 1i*randn(Nsrc,Mw,Nfreq))/sqrt(2); %TF of q(t) for Mw snapshots and Nfreq frequencies


for f=1:Nfreq

	%correlate sources
	Q(:,:,f)=Qcorr*Q(:,:,f);

	%choose a propagator
	%G_src_mic(:,:,f) = randn(Nmic,Nsrc)+1i*randn(Nmic,Nsrc); %random mixing matrix
    %G_src_mic(:,:,f) = GreenFreeFieldUniformFlow(r_src , r_mic , freqs(f) , -fftsign ,Mach); %monopole in uniform flow
    G_src_mic(:,:,f) = GreenFreeField(r_src, r_mic , freqs(f) , -fftsign); %monopole in free field
        
    %Calculate propagated signals
	P(:,:,f) = G_src_mic(:,:,f)*Q(:,:,f);
	
	%Generate Noise
	nrms = mean(abs(P(:,:,f))).*10.^((-SNR)./20); %rms value of the noise on each microphone.
	N(:,:,f) = nrms.*(randn(Nmic,Mw) + 1i*randn(Nmic,Mw))/sqrt(2);
	
	%extra noise
	nrms_extra=	mean(abs(P(:,:,f))).*10.^((-extra_SNR)./20);
	N(noisy_mic,:,f)= nrms_extra.*(randn(nb_noisy_mic,Mw) + 1i*randn(nb_noisy_mic,Mw))/sqrt(2);
	
	%Add noise to signals 	
	Y(:,:,f) = P(:,:,f) + N(:,:,f);
	
	%Calculate cross spectra 
	Sq(:,:,f) = Q(:,:,f) * Q(:,:,f)' /Mw;
	Sy(:,:,f) = Y(:,:,f) * Y(:,:,f)' /Mw; %signal + noise
	Sp(:,:,f) = P(:,:,f) * P(:,:,f)' /Mw; %mic signals
	Sn(:,:,f) = N(:,:,f) * N(:,:,f)' /Mw; %noise only
	%rank(Sy(:,:,f))
	%disp(['Rang de Sp : ' num2str(rank(Sp(:,:,f),0.005*max(eig(Sp(:,:,f))))) ]);	 
end
%%

% %%%--------------------------------------------------------------------------------------------
% %%% Display data
% %%%--------------------------------------------------------------------------------------------
% 
% figure
% plot3(r_mic(:,1),r_mic(:,2),r_mic(:,3),'ob')
% hold on
% plot3(x_src(:),y_src(:),z_src(:),'or')
% title('Configuration')
% legend('mic positions','src positions')
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')

%%

% %%%--------------------------------------------------------------------------------------------
% %%% Export data
% %%%--------------------------------------------------------------------------------------------
% Sp_filename='Sp_SNR0dB_Mw9300_49src_f3000.h5';
% 
% if exist(Sp_filename)
% 	d=input(['Delete ' Sp_filename ' (y/n) ? '],'s');
% 	
% 	if d=='y'
% 		delete(Sp_filename);
% 		h5create(Sp_filename,'/signal',size(Sp));
% 		h5write(Sp_filename,'/signal',Sp);
% 		
% 		h5create(Sp_filename,'/signal_plus_noise',size(Sy));
% 		h5write(Sp_filename,'/signal_plus_noise',Sy);
%         
%         h5create(Sp_filename,'/noise',size(Sn));
% 		h5write(Sp_filename,'/noise',Sn);
% 		
% 		h5create(Sp_filename,'/sources',size(Sq));
% 		h5write(Sp_filename,'/sources',Sq);
% 
% 		h5create(Sp_filename,'/frequencies',size(freqs));
% 		h5write(Sp_filename,'/frequencies',freqs);
%         
%         disp('Data saved');
% 	else 
% 		disp('No data saved');
% 	end
% else
% 		delete(Sp_filename);
% 		h5create(Sp_filename,'/signal',size(Sp));
% 		h5write(Sp_filename,'/signal',Sp);
% 		
% 		h5create(Sp_filename,'/signal_plus_noise',size(Sy));
% 		h5write(Sp_filename,'/signal_plus_noise',Sy);
%         
%         h5create(Sp_filename,'/noise',size(Sn));
% 		h5write(Sp_filename,'/noise',Sn);
% 
% 		h5create(Sp_filename,'/frequencies',size(freqs));
% 		h5write(Sp_filename,'/frequencies',freqs);
%         
%         disp('Data saved');
% 	
% end

