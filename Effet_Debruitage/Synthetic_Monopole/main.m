%Alice Dinsenmeyer, mars 2018
%Étude de l'influence du débruitage sur la quantification d'un monopole

clear all;
close all;

addpath('/home/adinsenmey/Bureau/these_2017/Effet_Debruitage/Generate_Spectra');

%%%-----------------------------------------------------------------------------
%%% Global constants
%%%-----------------------------------------------------------------------------
methods={'none'};%{'DRec','RPCA','PFA','DR'}; %denoising methods to apply
fftsign=+1;

%%%-----------------------------------------------------------------------------
%%% Experiment setup
%%%-----------------------------------------------------------------------------
freq=15000; %frequency
rho=0; %correlation between monopoles
SNR=10^10;%noise level added
Mw=10^4; %number of snapshots

%coordinates antenna
r_mic=load('spiral_array_coord.mat');
r_mic=r_mic.mic_info;
Nmic = size(r_mic,1);

%coordinates of the source
x_src=[0.4];%[0.4 1];
z_src=[0.8];%[0.5 1];
Nsrc=length(x_src); %number of monopole


%beamforming map coordinates
Nx=51;
Nz=51;
Ny=1;

y_grid=linspace(0,0,Ny);
x_grid=linspace(-1,1,Nx);
z_grid=linspace(-1,1,Nz);

[x_map,y_map, z_map]=meshgrid(x_grid,y_grid,z_grid);

r_map=[x_map(:) y_map(:) z_map(:)];

%Green function from map to antenna
G_map_mic = GreenFreeField(r_map, r_mic , freq, -fftsign); %monopole in free field



%%%-----------------------------------------------------------------------------
%%% Generate CSM
%%%-----------------------------------------------------------------------------
[Sq Sy Sp Sn ] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw,SNR , x_src , z_src); %no extra-noise

%%

for i=1:length(methods)

    %%%-----------------------------------------------------------------------------
    %%% Denoising
    %%%-----------------------------------------------------------------------------
    switch methods{i}
        case 'DRec',
            cvx_quiet('true'); cvx_precision('high');
            [Sy_denoised]=CSMRecHald(Sy);
            figure(66)
            hold on
            plot(real(diag(Sy_denoised)))
            
        case 'RPCA',
        
        case 'PFA',
            
        case 'DR'
            Sy_denoised=Sy;
            Sy_denoised(logical(eye(Nmic)))=0;
        
        case 'none' %no denoising
            Sy_denoised=Sy;
        
            
    end

    %%%-----------------------------------------------------------------------------
    %%% Beamforming
    %%%-----------------------------------------------------------------------------
    dist_norm = (sqrt(sum(abs((r_map-mean(r_mic))).^2,2)))';%distance map -> centre de l'antenne
    [BF_map_diag] = beamforming(diag(diag(Sy_denoised)) , G_map_mic);
    %BF_map_diag = BF_map_diag./dist_norm; %normalisation par la distance au centre de l'antenne
    BF_map_diag=reshape(BF_map_diag,Nx,Nz);

    [BF_map_extra] = beamforming(Sy_denoised-diag(diag(Sy_denoised)) , G_map_mic);
    %BF_map_extra=BF_map_extra./dist_norm;
    BF_map_extra=reshape(BF_map_extra,Nx,Nz);
    
    BF_map = BF_map_diag + BF_map_extra;
    
    %set negative values to zero for log
    BF_map_diag(BF_map_diag<0)=0;
    BF_map_extra(BF_map_extra<0)=0;
    BF_map(BF_map<0)=0;
    
    %%%-----------------------------------------------------------------------------
    %%% Display results
    %%%-----------------------------------------------------------------------------
    
    figure
    suptitle(['Beamforming map, denoising method: ' methods{i}])
    subplot(1,3,1)
    imagesc(x_grid,z_grid,(real(BF_map_diag)));
    %caxis([-30 5])
    set(gca,'Ydir','Normal')

    
    subplot(1,3,2)
    imagesc(x_grid,z_grid,(real(BF_map_extra)));
    %caxis([-30 5])
	set(gca,'Ydir','Normal')

    
    subplot(1,3,3)
    imagesc(x_grid,z_grid,(real(BF_map)));
    %caxis([-30 5])
    set(gca,'Ydir','Normal')

        
    %%%-----------------------------------------------------------------------------
    %%% Save Data
    %%%-----------------------------------------------------------------------------
    
end
