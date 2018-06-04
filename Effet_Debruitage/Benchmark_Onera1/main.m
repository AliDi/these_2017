% Alice Dinsenmeyer mars-avril 2018
% Etude de l'effet du débruitage sur le niveau de source retrouvé par des
% méthodes d'imagerie : 
%-beamforming
%
% Etude appliquée aux données du benchmark Onera1 (open section)
clear all; 
%close all;

addpath('/home/adinsenmey/Bureau/these_2017/Effet_Debruitage/Generate_Spectra');
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

%%%------------------------------------------------------------------------
%%% Chargement des données
%%%------------------------------------------------------------------------

data_path = '/media/adinsenmey/ECHANGE/BENCHMARK/experimental/ONERA1/Fourniture/';
case_number=6; %selon la notation de la doc du benchmark.

switch case_number
    case 1
        file_name = 'BenchmarkOnera_SingleSourceWithoutFlow_pt0032_CSMCorr.h5';     %1 source, M=0
        Nsrc=1;
    case 5
        file_name = 'BenchmarkOnera_SingleSourceWithFlow_pt0020_CSMCorr.h5';        %1 source, M=0.13
        Nsrc=1;
    case 3
        file_name = 'BenchmarkOnera_TwoSourceCorrPPWithoutFlow_pt0036_CSMCorr.h5';  %2 sources corrélées (en phase), M=0
        Nsrc=2;
    case 7
        file_name = 'BenchmarkOnera_TwoSourceCorrPPWithFlow_pt0021_CSMCorr.h5';     %2 sources corrélées (en phase), M=0.13
        Nsrc=2;
    case 2
        file_name = 'BenchmarkOnera_TwoSourceUncorrWithoutFlow_pt0038_CSMCorr.h5';  %2 sources décorrélées, M=0
        Nsrc=2;
    case 6
        file_name = 'BenchmarkOnera_TwoSourceUncorrWithFlow_pt0024_CSMCorr.h5';     %2 sources décorrélées, M=0.13
        Nsrc=2;

end

%%

%Misc
fftsign=single(h5readatt([data_path file_name], '/CsmData','fftSign'));
freqs=h5read([data_path file_name],'/CsmData/binCenterFrequenciesHz');
num_freq=length(freqs);
Mach = h5read([data_path file_name],'/MeasurementData/machNumber');
c = h5read([data_path file_name],'/MeasurementData/speedOfSoundMPerS');
Mw=600; %approximate number of snapshots

% Mic Info
num_mic=h5readatt([data_path file_name],'/MetaData/ArrayAttributes','microphoneCount');
mic_info=h5read([data_path file_name],'/MetaData/ArrayAttributes/microphonePositionsM');
coor_ref=h5readatt([data_path file_name],'/MetaData/TestAttributes','coordinateReference');

%CSM
cpimag=h5read([data_path file_name],'/CsmData/csmImaginary');
cpreal=h5read([data_path file_name],'/CsmData/csmReal');
Sy_non_reshape=cpreal + i*cpimag;
for f=1:num_freq
    Sy(:,:,f)= Sy_non_reshape(f,:,:);    
end
Sy=double(Sy);


%Eclusion du 71e et 115e mic
Sy=Sy([1:70 72:end-1],[1:70 72:end-1],:);
mic_info = mic_info([1:70 72:end-1],:);
num_mic=size(mic_info,1);


%freq_study =round([2500 5000 10000 20000]./(freqs(2)-freqs(1)))+1;
freq_study =round([10000]./(freqs(2)-freqs(1)))+1;


%Setup

r_L2 = [2.25 -0.11 -0.03]; %position de la source L2
r_L4 = [2.41 -0.11 -0.03]; %position de la source L4

n_SL = [0 1 0]; %vecteur normal à l'interface plane ecoulement/salle
pos_SL = [0 -1.6 0]; %point de repérage de la position de la couche limite



%%

%%%------------------------------------------------------------------------
%%% DEFINE GRID FOR SOURCE MAP AND TRANSFERT MATRIX
%%%------------------------------------------------------------------------
Nx=80;
Ny=1;
Nz=30;
Nmap = Nx * Ny * Nz; %number of source point on the constructed map

grid_bounds=h5read([data_path file_name],'/MetaData/TestAttributes/domainBoundsM');


x_grid = linspace(grid_bounds(1,1),grid_bounds(2,1),Nx)';
y_grid = linspace(grid_bounds(1,2),grid_bounds(2,2),Ny)';
z_grid = linspace(grid_bounds(2,3)-0.1 , grid_bounds(1,3)+0.1 , Nz)';

dx_grid = x_grid(2)-x_grid(1);
dz_grid = z_grid(2)-z_grid(1);


[Xmap,Zmap, Ymap] = meshgrid(x_grid,z_grid,y_grid);
%Xmap=Xmap'; Ymap=Ymap'; Zmap=Zmap';

r_src_map = [Xmap(:) Ymap(:) Zmap(:)]; % vector of estimate source point positions (map)

%%% Matrix of acoustic transfert from map points to mic position (for the BF pb)
i=1;
for f=freq_study
    disp(['Gmapmic : Frequency : ' num2str(f) '/' num2str(num_freq)]);
    G_map_mic(:,:,i) = GreenFreeFieldInfShearLayer(r_src_map , mic_info , freqs(f) , fftsign ,Mach , c , n_SL , pos_SL );
    %G_map_mic(:,:,i) = GreenFreeFieldUniformFlow(r_src_map , mic_info , freqs(f) , fftsign ,Mach,c);
    %G_map_mic(:,:,i) = GreenFreeField(r_src_map, mic_info , freqs(f) , fftsign , c);
    i=i+1;
end

%%
%%%-----------------------------------------------------------------------------
%%% Denoising
%%%-----------------------------------------------------------------------------
%Parameters for PFA
K_est=5;
Nrun=1000;
opt.noise='hetero';

%Parameters for NNLS
options=optimset('TolX',1e-15,'MaxIter',2000); %set termination tolerance for Sqq

methods={'none','DR','PFA','RPCA','DRec','DR_Quentin'};%,'DR','PFA'};

i=1;
for f=freq_study

    for j=1:length(methods)
    
        flag_bf_quentin=0; %no need to apply Quentin's BF

		switch methods{j}
		    case 'DRec',
		        cvx_quiet('true'); cvx_precision('high');
		        [Sy_denoised]=CSMRecHald(Sy(:,:,f));
		        
                %{
                figure(66)
		        hold on
		        plot(real(diag(Sy_denoised)))
                %}
		        
		    case 'RPCA',
                [Sy_denoised] = proximal_gradient_rpca(Sy(:,:,f) , 1/sqrt(num_mic), 300,1e-7,-1,-1,-1,-1);
		    
		    case 'PFA',
                a=struct();
                b=struct();
                Ini=struct();
    
                noise = real(mean(diag(Sy(:,:,f))));
	            alpha2_mean = abs(real(sort(eig(Sy(:,:,f))/max(eig(Sy(:,:,f))),'descend')));
	   
                %%% hyper-parametres

                %alpha
                a.alpha = 1./alpha2_mean(1:K_est);      % hyper-hyper-paramètres sur alpha
	
                %beta
                [a.beta2,b.beta2] = Convert_InvGamma(mean(noise.^2),10*mean(noise.^2));     % hyper-hyper-paramètres sur beta2
	
                %gamma
                gamma_mean = (real(trace(Sy(:,:,f)))/num_mic )/mean(alpha2_mean); %Alice : pas de bruit retiré
                [a.gamma2,b.gamma2] = Convert_InvGamma(gamma_mean,10*gamma_mean);   % hyper-hyper-paramètres sur gamma2    
                
                %initialisation
                Ini.alpha(1,:) = 1./a.alpha(:)';
                Ini.beta2(1,:) = b.beta2/a.beta2*ones(num_mic,1);
                Ini.gamma2(1) = b.gamma2/a.gamma2;
                Ini.Lambda(1,:,:) = (randn(num_mic,K_est) + 1i*randn(num_mic,K_est))/sqrt(2);

                %opt.gamma2=1.1;    
                [Sc,Lambda,alpha,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3(Sy(:,:,f),K_est,a,b,Mw,Nrun,opt,Ini);
                for jj=1:Nrun
                    tmp_Sy_denoised(:,:,jj)=squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))';
                end
                Sy_denoised=mean(tmp_Sy_denoised(:,:,500:end),3);
		        
		    case 'DR'
		        Sy_denoised=Sy(:,:,f);
		        Sy_denoised(logical(eye(num_mic)))=0;
		    
		    case 'none' %no denoising
		        Sy_denoised=Sy(:,:,f);
		        
		    case 'DR_Quentin'
                flag_bf_quentin=1; %apply Quentin's BF
		    
		        
		end

	%%%------------------------------------------------------------------------
	%%% Beamforming
	%%%------------------------------------------------------------------------



    disp(['BF : Frequency : ' num2str(f) '/' num2str(num_freq)]);
    a=zeros(Nx*Ny*Nz,1);
    b=zeros(Nx*Ny*Nz,num_mic);
    
    if flag_bf_quentin
        [a b] = beamforming_DR_Quentin(Sy(:,:,f)-diag(diag(Sy(:,:,f))) , G_map_mic(:,:,i));
    else
    	[a b] = beamforming(Sy_denoised , G_map_mic(:,:,i) );
    end
     
    %[a , ~ , count ] = clean_sc(1 , Sy_denoised , G_map_mic(:,:,i) , 300 , 'on');


    
    %%%-----------------------------------------------------------------------------
%%  % Deconvolution
	%%%-----------------------------------------------------------------------------
	[a] = nnls(double(real(a)) , double(G_map_mic(:,:,i)) , double(b), options);
    
     %%%-----------------------------------------------------------------------------
%%  % Save data
	%%%-----------------------------------------------------------------------------
    amp_src = sort(real(a),'descend');
    amp_src = amp_src(1:Nsrc)';
    
    disp('=================================')
    disp([ num2str(methods{j}) ' -> Amp (dB) = ' num2str(10*log10(amp_src(:)')) ])
    disp('=================================')

    
	%%%-----------------------------------------------------------------------------
	%%% Display Beamforming maps
	%%%-----------------------------------------------------------------------------
    %%%reshape map in2D
    Sqq_BF(:,:,i)=reshape(real(a),Nz,Nx);
    
    %%%display results
    A=real(Sqq_BF(:,:,i));
    A(A<0)=0;
    A=20*log10(sqrt(A));
    
	figure
	imagesc(x_grid,z_grid,A);
	colorbar
	xlabel('x (m)');
	ylabel('y (m)');
	title(['Beamforming, f=' num2str(freqs(f)) ', d\''ebruitage :' num2str(methods{j}) ]) 
    hold on
    plot(r_L2(1),r_L2(3),'ro')
    plot(r_L4(1),r_L4(3),'ro') 
    caxis([max(max(A))-15 max(max(A))])
    
    end
    i=i+1;


end


