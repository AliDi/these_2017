%Beamforming on 1D configuration (ie mics aligned) : 
%M mics
%nb_F frequencies
%N sources


clear all ;
%close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')


%%%------------------------------------------------------------------------
%%% CONSTANT PARAMETERS
%%%------------------------------------------------------------------------

%%% Frequencies
global F;
F=1:100:20e3; %observed frequencies
nb_F = length(F);

%%% Microphones position
array_length = 0.7; %in meters
M = 64; %number of microphones
z_mic = 0 * ones(1,M);
x_mic = linspace(0,0.7,M);%x_mic = sort([0.5462    0.2728    0.1692    0.2827    0.0675    0.0924    0.6594    0.6693    0.4026    0.0418    0.1643    0.2472  0.5748    0.0108    0.0301    0.1183    0.4544    0.5122    0.4534    0.3156    0.3829    0.2074    0.5213    0.1323   0.4807    0.1285    0.2579    0.4379    0.5462    0.0568    0.6506    0.5430 ]) ;%sort(array_length*rand(1, M));
y_mic = zeros(1,M);
r_mic = [x_mic' y_mic' z_mic'];

%%% Real sources position
N=14; %number of sources
z_src = 0.7 * ones(1,N);
x_src =logspace(log10(0.10),log10(0.60),N); 
y_src = zeros(1,N);
r_src = [x_src' y_src' z_src'];

figure(1)
plot(x_mic,z_mic,'ob')
hold on
plot(x_src,z_src,'or')
title('Configuration')

%%% Type of source : correlated (correlated = 'correlated') or uncorrelated (correlated = 'uncorrelated')
correlated ='correlated'

%%Real source strength coefficients
q = 1*ones(N , 1 , nb_F); 

if strcmpi(correlated,'correlated')
	Sqq=InterspectraMatrix(q);
	
else
	for f=1:nb_F
		Sqq(:,:,f)=diag(q(:,:,f)); %source inter spectra is null
	end
end



%%%------------------------------------------------------------------------
%%% DIRECT PB : generates mic pressures
%%%------------------------------------------------------------------------

%%% Define grid for source map
Nx=1000;
Ny=1;
Nz=1;
Nmap = Nx * Ny * Nz; %number of source point on the constructed map

x_grid = linspace(-0.25,0.85,Nx);
y_grid = linspace(0,0,Ny);
z_grid = linspace(0.7 , 0.7 , Nz);

[Xmap,Ymap, Zmap] = meshgrid(x_grid,y_grid,z_grid);
Xmap=Xmap';
Ymap=Ymap';
Zmap=Zmap';

r_src_map = [Xmap(:) Ymap(:) Zmap(:)]; % vector of estimate source point positions (map)

%%% Matrix of acoustic transfert from real source position to mic position (for the direct pb)
G_src_mic = GreenFreeField(r_src , r_mic);




for f=1:nb_F

	
	%%% Direct problem : generate mic pressures
	Spp(:,:,f) = G_src_mic(:,:,f)*Sqq(:,:,f)*G_src_mic(:,:,f)';
end



%%%------------------------------------------------------------------------
%%% BEAMFORMING
%%%------------------------------------------------------------------------

%%% Matrix of acoustic transfert from map points to mic position (for the BF pb)
G_map_mic = GreenFreeField(r_src_map , r_mic);

for f=1:nb_F
	
	%%% Calculate weights/source scaling
	for m=1:M
		Lvect(m)=norm(G_map_mic(:,m,f));
	end	
	L=diag(Lvect);
	
	%%% Calculate sterring vectors
    for n=1:Nmap
		W(n,:,f)=G_map_mic(:,n,f)'*L;
        %W(n,:,f)=G_map_mic(:,n,f)'.*(G_map_mic(:,n,f)'*G_map_mic(:,n,f))^(-1);
    end
	
	%%% Estimates sources  
	Sqq_est(:,:,f)=W(:,:,f)*Spp(:,:,f)*W(:,:,f)';
	Sqq_diag(:,f)=diag(Sqq_est(:,:,f)); %save auto spectra
	
end


%%%-----------------------------------------------------------------------
%%% Display results
%%%-----------------------------------------------------------------------
% sources on a line (1D)
figure(4)
imagesc(F,x_grid,20*log10(abs(Sqq_diag)));
colorbar
xlabel('Fr\''equence en Hz');
ylabel('Position en m');
title(['Beamforming, ' correlated ' sources'])

print(['beamforming_' correlated],'-dsvg')
print(['beamforming_' correlated],'-dpng')






