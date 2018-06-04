%Alice Dinsenmeyer, mars 2018
%Étude de l'influence du débruitage sur la quantification d'un monopole

clear all;
%close all;

addpath('/home/adinsenmey/Bureau/these_2017/Effet_Debruitage/Generate_Spectra');
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')

%%%-----------------------------------------------------------------------------
%%% Global constants
%%%-----------------------------------------------------------------------------
methods={'none'}%{'RPCA','DR','DR_Quentin'};%{'DRec','RPCA','DR','DR_Quentin','none'};%{'none','DRec','RPCA','PFA','DR'}; %denoising methods to apply
fftsign=+1;

%%%-----------------------------------------------------------------------------
%%% Experiment setup
%%%-----------------------------------------------------------------------------
freq=4000; %frequency
rho=0; %correlation between monopoles
SNR=-20;%noise level added
Mw=10^4; %number of snapshots
c=340;

%coordinates antenna
r_mic=load('spiral_array_coord.mat');
r_mic=r_mic.mic_info;


%antenne rectangle regulière
% [tmp1 tmp2]=meshgrid(linspace(-1,1,10), linspace(-1,1,9));
% r_mic(:,1)=tmp1(:);
% r_mic(:,3)=tmp2(:);

Nmic = size(r_mic,1);

%sources
y_src=0;
Nsrc_study=unique(round(logspace(log10(1),log10(10),10)));




%beamforming map coordinates
Nx=25;
Nz=61;
Ny=1;

x_grid_min=-0.5;
x_grid_max=0.5;
z_grid_min=-1.2;
z_grid_max=1.2;
y_grid=linspace(0,0,Ny);
x_grid=linspace(x_grid_min,x_grid_max,Nx);
z_grid=linspace(z_grid_min,z_grid_max,Nz);

dx_grid = abs(x_grid_min-x_grid_max)/(Nx-1);
dz_grid = abs(z_grid_min-z_grid_max)/(Nz-1);


[x_map,y_map, z_map]=meshgrid(x_grid,y_grid,z_grid);

r_map=[x_map(:) y_map(:) z_map(:)];



%Green function from map to antenna
G_map_mic = GreenFreeField(r_map, r_mic , freq, -fftsign,c); %monopole in free field

%Delete former data files and create new file
for i=1:length(methods)
    file_name=['out_data/' methods{i} '_Nsrc.h5'];
    delete(file_name);
    h5create(file_name,'/Nsrc', [length(Nsrc_study) 1]);
    h5create(file_name,'/Amplitude_Source',[length(Nsrc_study) max(Nsrc_study)]);
end

%Parameters for PFA
Nrun=1000;
opt.noise='hetero';

%Parameters for NNLS
options=optimset('TolFun',1e-3,'TolX',1e-3,'MaxIter',2000); %set termination tolerance for Sqq
    
%%

for j=1:length(Nsrc_study)
    disp(['Nsrc = ' num2str(Nsrc_study(j))])
    
    %coordinates of the source
	%one line with rotation theta
	Nsrc=Nsrc_study(j);
	theta=0.0175;
	R=0.9;
	z_src = linspace(-R,R,Nsrc);
    x_src=zeros(1,Nsrc);
    
	x_src_rot =x_src*cos(theta) - z_src*sin(theta) +0.2;
	z_src_rot=x_src*sin(theta)+z_src*cos(theta) ;    
    x_src=round(x_src_rot/dx_grid)*dx_grid;
    z_src=round(z_src_rot/dz_grid)*dz_grid;
    
    y_src=0;
	[x_src y_src]=meshgrid(x_src,y_src);
	r_src = [x_src(:) y_src(:) z_src(:)]; %coordinates of real sources
		
    
    %Check if source is on a grid point
	ind_x_src=(x_src-x_grid_min)/dx_grid + 1;%indice de la grille correspondant à la position de la source
	ind_z_src=(z_src-z_grid_min)/dz_grid + 1;
	if sum((double(int32(ind_x_src))-ind_x_src)> 1e-5) ||sum((double(int32(ind_z_src))-ind_z_src)> 1e-5)
		error('Une des sources est pas sur un pas de grille')
    end    
    
    %%%-----------------------------------------------------------------------------
    %%% Generate CSM
    %%%-----------------------------------------------------------------------------
    [Sq Sy Sp Sn ] = generate_Spp_signal(freq, Nsrc , rho , SNR , Mw,SNR , x_src , z_src,r_mic); %no extra-noise
    %d_ref=real(diag(Sp));

    
    for i=1:length(methods)
        disp(['Denoising method : ' methods{i}]);
        
        flag_bf_quentin=0; %no need to apply Quentin's BF


        %%%-----------------------------------------------------------------------------
        %%% Denoising
        %%%-----------------------------------------------------------------------------
        switch methods{i}
            case 'DRec',
                cvx_quiet('true'); cvx_precision('high');
                [Sy_denoised]=CSMRecHald(Sy);

                %{
                figure(66)
                hold on
                plot(real(diag(Sy_denoised)))
                %}

            case 'RPCA',
                [Sy_denoised] = proximal_gradient_rpca(Sy , 1/sqrt(Nmic), 300,1e-7,-1,-1,-1,-1);

            case 'PFA',
                K_est=Nsrc+5;
                if K_est>Nmic
                    K_est=Nmic;
                end
                a=struct();
                b=struct();
                Ini=struct();

                noise = real(mean(diag(Sy)));
                alpha2_mean = abs(real(sort(eig(Sy)/max(eig(Sy)),'descend')));

                %%% hyper-parametres

                %alpha
                a.alpha = 1./alpha2_mean(1:K_est);      % hyper-hyper-paramètres sur alpha

                %beta
                [a.beta2,b.beta2] = Convert_InvGamma(mean(noise.^2),10*mean(noise.^2));     % hyper-hyper-paramètres sur beta2

                %gamma
                gamma_mean = (real(trace(Sy))/Nmic )/mean(alpha2_mean); %Alice : pas de bruit retiré
                [a.gamma2,b.gamma2] = Convert_InvGamma(gamma_mean,10*gamma_mean);   % hyper-hyper-paramètres sur gamma2    

                %initialisation
                Ini.alpha(1,:) = 1./a.alpha(:)';
                Ini.beta2(1,:) = b.beta2/a.beta2*ones(Nmic,1);
                Ini.gamma2(1) = b.gamma2/a.gamma2;
                Ini.Lambda(1,:,:) = (randn(Nmic,K_est) + 1i*randn(Nmic,K_est))/sqrt(2);

                %opt.gamma2=1.1;    
                [Sc,Lambda,alpha,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3(Sy,K_est,a,b,Mw,Nrun,opt,Ini);
                tmp_Sy_denoised=zeros(Nmic,Nmic,Nrun);
                for jj=1:Nrun
                    tmp_Sy_denoised(:,:,jj)=squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))';
                    %tmp_d_mcmc(:,jj)=real(diag(squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))'));
                    %tmp_err_mcmc(jj,1) = norm( d_ref - tmp_d_mcmc(:,jj)) /  norm( d_ref);        
                end
                %err_mcmc(j,:)=tmp_err_mcmc;
                Sy_denoised=mean(tmp_Sy_denoised(:,:,500:end),3);

            case 'DR'
                Sy_denoised=Sy;
                Sy_denoised(logical(eye(Nmic)))=0;

            case 'none' %no denoising
                Sy_denoised=Sy;
             
            case 'DR_Quentin'
                flag_bf_quentin=1; %apply Quentin's BF
        end

		%%%-----------------------------------------------------------------------------
		%%% Beamforming
		%%%-----------------------------------------------------------------------------
		 if flag_bf_quentin
            [BF_map , W] = beamforming_DR_Quentin(Sy-diag(diag(Sy)) , G_map_mic);
            BF_map=reshape(BF_map,Nx,Nz);
            
        else %apply standart BF 
        
        %{           
            [BF_map_diag] = beamforming(diag(diag(Sy_denoised)) , G_map_mic);
            BF_map_diag=reshape(BF_map_diag,Nx,Nz);
            
            [BF_map_extra] = beamforming(Sy_denoised-diag(diag(Sy_denoised)) , G_map_mic);
            BF_map_extra=reshape(BF_map_extra,Nx,Nz);

            BF_map = BF_map_diag + BF_map_extra;
       %}
            
            [BF_map ,W ] = beamforming(Sy_denoised , G_map_mic);
            
            
            %---BF_map_diag = BF_map_diag./dist_norm; %normalisation par la distance au centre de l'antenne
            %---dist_norm = (sqrt(sum(abs((r_map-mean(r_mic))).^2,2)))';%distance map -> centre de l'antenne
            %---BF_map_extra=BF_map_extra./dist_norm;
        end
        
        %%%-----------------------------------------------------------------------------
		%%% Deconvolution
		%%%-----------------------------------------------------------------------------
		[BF_map] = nnls(real(BF_map) , G_map_mic , W , options);
        BF_map=reshape(BF_map,Nx,Nz);

        imagesc(real(BF_map))
        pause(0.5)
%%		
%{	
		
		%%%-----------------------------------------------------------------------------
		%%% Display Beamforming maps
		%%%-----------------------------------------------------------------------------
		
        %set negative values to zero for log
		%BF_map_diag(BF_map_diag<0)=0;
		%BF_map_extra(BF_map_extra<0)=0;
		%BF_map(BF_map<0)=0;
        
		figure
		suptitle(['Beamforming map, denoising method: ' methods{i}])
		subplot(1,3,1)
		imagesc(x_grid,z_grid,(real(BF_map_diag')));
		%caxis([-30 5])
        xlabel('x (m)')
        ylabel('z (m)')
        title('BF(diag($\bm{S_{pp}}$))')
        hold on
        plot(x_src,z_src,'*r')
        colorbar
        set(gca,'Ydir','Normal')
        %caxis([real(min(min(BF_map_extra))) real(max(max(BF_map_extra)))])
        
		
		subplot(1,3,2)
		imagesc(x_grid,z_grid,(real(BF_map_extra')));
		%caxis([-30 5])
		xlabel('x (m)')
        ylabel('z (m)')
        title('BF($\bm{S_{pp}}$-diag($\bm{S_{pp}}$))')
        hold on
        plot(x_src,z_src,'*r')
        colorbar
        caxis([real(min(min(BF_map_extra))) real(max(max(BF_map_extra)))])
        set(gca,'Ydir','Normal')

		
		subplot(1,3,3)
		imagesc(x_grid,z_grid,(real(BF_map')));
		%caxis([-30 5])
        xlabel('x (m)')
        ylabel('z (m)')
        title('BF($\bm{S_{pp}}$)')
        hold on
        plot(x_src,z_src,'*r')
        colorbar
        caxis([real(min(min(BF_map_extra))) real(max(max(BF_map_extra)))])
		set(gca,'Ydir','Normal')

		
		
%}	
%%
		    
		%%%-----------------------------------------------------------------------------
		%%% Save Data
		%%%-----------------------------------------------------------------------------
        
        %récupère les amplitudes à l'emplacement des vraies sources
        %{
 		for s=1:Nsrc
             amp_src(s) = BF_map(int32(ind_x_src(s)),int32(ind_z_src(s)));
        end
        %}
        
        %récupère les Nsrc plus grandes amplitudes de source
        amp_src = sort(BF_map(:),'descend');
        amp_src = amp_src(1:Nsrc)';
                
        h5write(['out_data/' methods{i} '_Nsrc.h5'],'/Nsrc',Nsrc,[j 1],[1 1]);
        h5write(['out_data/' methods{i} '_Nsrc.h5'],'/Amplitude_Source',amp_src,[j 1],[1 Nsrc]);
       		
	end
end

%%
%%%-----------------------------------------------------------------------------
%%% Display results
%%%-----------------------------------------------------------------------------
figure


for i=1:length(methods)
        amp=h5read(['out_data/' methods{i} '_Nsrc.h5'],'/Amplitude_Source');
        Nsrc=h5read(['out_data/' methods{i} '_Nsrc.h5'],'/Nsrc');
        semilogx(Nsrc_study,10*log10(sum(amp(:,:),2)./Nsrc),'o-','DisplayName',methods{i})
        hold on
end
semilogx(Nsrc_study,zeros(1,length(Nsrc_study)),':')
xlabel('Nsrc')
ylabel('$10 \log_{10}(\tilde{S_{qq_0}}/S_{qq_0})$')
    
legend('show');


