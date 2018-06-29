%Alice Dinsenmeyer, mars 2018
%Étude de l'influence du débruitage sur la quantification d'un monopole

clear all;
%close all;

addpath('/home/adinsenmey/Bureau/these_2017/Effet_Debruitage/Generate_Spectra');
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/imagerie_jerome')

%%%-----------------------------------------------------------------------------
%%% Global constants
%%%-----------------------------------------------------------------------------
methods={'none'};%{'none','DRec','RPCA','PFA','DR'}; %denoising methods to apply
BF_flag=1; %display Beamforming results or not
fftsign=+1;

%%%-----------------------------------------------------------------------------
%%% Experiment setup
%%%-----------------------------------------------------------------------------
freq=1200; %frequency
rho=0.999; %correlation between monopoles
SNR=0;%noise level added
Mw=10^4; %number of snapshots
c=340;
lambda=c/freq;
k=2*pi/lambda;

%coordinates antenna
r_mic=load('spiral_array_coord.mat');
r_mic=r_mic.mic_info(1:93,:);
r_mic(:,2)=r_mic(:,2);

%antenne rectangle regulière
% [tmp1 tmp2]=meshgrid(linspace(-1,1,10), linspace(-1,1,9));
% r_mic(:,1)=tmp1(:);
% r_mic(:,3)=tmp2(:);

Nmic = size(r_mic,1);

%beamforming map coordinates
Ny=1;

dx_grid=lambda/7;
dz_grid=dx_grid;

x_grid_min=-8*dx_grid;
x_grid_max=8*dx_grid;
z_grid_min=-8*dz_grid;
z_grid_max=8*dz_grid;

x_grid=x_grid_min:dx_grid:x_grid_max;
Nx=length(x_grid);

z_grid=z_grid_min:dz_grid:z_grid_max;
Nz=length(z_grid);

Nmap=Nx*Nz;

y_grid=linspace(0,0,Ny);
%x_grid=linspace(x_grid_min,x_grid_max,Nx);
%z_grid=linspace(z_grid_min,z_grid_max,Nz);

%dx_grid = abs(x_grid_min-x_grid_max)/(Nx-1);
%dz_grid = abs(z_grid_min-z_grid_max)/(Nz-1);


[x_map,y_map, z_map]=meshgrid(x_grid,y_grid,z_grid);

r_map=[x_map(:) y_map(:) z_map(:)];

%coordinates of the source
%one line with rotation theta
Nsrc=1;
theta=0;%.0175;
R=0.07;
z_src = linspace(-R,R,Nsrc);
x_src=zeros(1,Nsrc);

x_src_rot =x_src*cos(theta) - z_src*sin(theta);
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
    error('La source n est pas sur un pas de grille')
end

%Green function from map to antenna
G_map_mic = GreenFreeField(r_map, r_mic , freq, -fftsign, c); %monopole in free field
%%

%Delete former data files and create new file
delete 'out_data/*';


for i=1:length(methods)
    file_name.map{i}=['out_data/MAP' methods{i} '_SNR.h5'];
    file_name.bf{i}=['out_data/BF' methods{i} '_SNR.h5'];
    
    %EM_map
    h5create(file_name.map{i},'/SNR', [length(SNR) 1]);
    h5create(file_name.map{i},'/Amplitude_Source',[length(SNR) Nx*Nz Nx*Nz]);
    h5create(file_name.map{i},'/Puissance',[1 length(SNR)]);

    
    %Beamforming
    h5create(file_name.bf{i},'/SNR', [length(SNR) 1]);
    h5create(file_name.bf{i},'/Amplitude_Source',[length(SNR) Nx Nz]);
    h5create(file_name.bf{i},'/Puissance',[1 length(SNR)]);
end

%Parameters for PFA denoising
K_est=Nsrc+2;%+5;
Nrun_den=100;
opt.noise='hetero';


%Parameters for MCMC FA
Nrun = 100;
opt = 'hetero';
   
%%

for j=1:length(SNR) %vary SNR
    disp(['SNR : ' num2str(SNR(j))]);

	%%%-----------------------------------------------------------------------------
	%%% Generate CSM
	%%%-----------------------------------------------------------------------------
	[Sq Sy Sp Sn ] = generate_Spp_signal(freq, Nsrc , rho , SNR(j) , Mw,SNR(j) , x_src , z_src,r_mic); %no extra-noise
    Sy=(Sy+Sy')/2;
	%%

    for i=1:length(methods)
		disp(['Denoising method : ' methods{i}]);
        


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
                [Sc,Lambda,alpha,beta2,gamma2] = MCMC_AnaFac_Quad_Sparse3(Sy,K_est,a,b,Mw,Nrun_den,opt,Ini);
                for jj=1:Nrun_den
                    tmp_Sy_denoised(:,:,jj)=squeeze(Lambda(jj,:,:))*diag(alpha(jj,:))*squeeze(Sc(jj,:,:))*diag(alpha(jj,:))*squeeze(Lambda(jj,:,:))';
                end
                Sy_denoised=mean(tmp_Sy_denoised(:,:,500:end),3);
		        
		    case 'DR'
		        Sy_denoised=Sy;
		        Sy_denoised(logical(eye(Nmic)))=0;
		    
		    case 'none' %no denoising
		        Sy_denoised=Sy;

		    
		        
        end
        Sy_denoised=(Sy_denoised+Sy_denoised')/2;
        %[Sy_denoised,n,errec, diag_matout]=recdiag(Sy_denoised,1,1000,1e-9,30);

        err_denoising = norm( real(diag(Sp)) - real(diag(Sy_denoised))) /  norm( real(diag(Sp)));  
        disp(['Erreur de débruitage : ' num2str(10*log10(err_denoising)) ' dB'])
%%
		%%%-----------------------------------------------------------------------------
		%%% Beamforming
		%%%-----------------------------------------------------------------------------
		%apply standart BF 
        if BF_flag==1   
        	[BF_map ,W ] = beamforming(Sy_denoised , G_map_mic);
            BF_map=reshape(BF_map,Nx,Nz);
            
            %---BF_map_diag = BF_map_diag./dist_norm; %normalisation par la distance au centre de l'antenne
            %---dist_norm = (sqrt(sum(abs((r_map-mean(r_mic))).^2,2)))';%distance map -> centre de l'antenne
            %---BF_map_extra=BF_map_extra./dist_norm;
		
			
			% Display Beamforming maps		
            
		    %set negative values to zero for log
			%BF_map_diag(BF_map_diag<0)=0;
			%BF_map_extra(BF_map_extra<0)=0;
			%BF_map(BF_map<0)=0;
		    
			figSrc=figure;
            subplot(1,2,1)
            A=(real(BF_map'));
            A(A<0)=0;
			imagesc(x_grid,z_grid,10*log10(A));
			%caxis([-30 5])
		    xlabel('x (m)')
		    ylabel('z (m)')
            title(['BF, denoising = ' methods{i}])
		    hold on
		    plot(x_src,z_src,'*r')
		    colorbar
            
            borne=max(max(10*log10(A)));
		    caxis([borne-20 borne])
			set(gca,'Ydir','Normal')
            
        end
%%
		%%%-----------------------------------------------------------------------------
		%%% MCMC FA
		%%%-----------------------------------------------------------------------------
        %%
        [u s v]=svd(diag(BF_map(:)));
        s=diag(s);
        
         noise = real(mean(diag(Sy_denoised)));
        %%% hyper-parametres

        %beta
        [a.beta2,b.beta2] = Convert_InvGamma(mean(noise),10*mean(noise));     % hyper-hyper-paramètres sur beta2

        %gamma
        alpha2_mean = abs(real(sort(eig(Sy_denoised),'descend'))); 
        [a.alpha2,b.alpha2] = Convert_InvGamma(mean(alpha2_mean),10*mean(alpha2_mean));   % hyper-hyper-paramètres sur alpha2    

        %initialisation
        Ini.alpha2(1,:) =alpha2_mean(1:K_est);%alpha2_mean(1:K_est);
        Ini.beta2(1,:) = b.beta2/a.beta2*ones(Nmic,1);
        
        Ini.Lambda(1,:,:) =(ones(Nmap,K_est)+1i*ones(Nmap,K_est));%
        
		[Sc_est,Lambda_est,alpha2_est,beta2_est] = MCMC_LatVar_QuadH(Sy_denoised,G_map_mic,K_est,a,b,Mw,Nrun,opt,Ini);
        %%
        for n=1:Nrun
        	tmp(:,:,n)=squeeze(Lambda_est(n,:,:))*squeeze(Sc_est(n,:,:))*squeeze(Lambda_est(n,:,:))';
            tmp2(:,n)=diag(tmp(:,:,n));
            %tmp(:,n)=diag(squeeze(Lambda_est(n,:,:))*diag(alpha2_est(n,:))*squeeze(Lambda_est(n,:,:))');
        end
        Sq_est_corr=mean(tmp(:,:,round(Nrun/2):end),3);        
        Sq_est=mean(tmp2(:,round(Nrun/2):end),2);
        
        %MAP.Sc(MAP.Sc<0)=0;
        Sq_est=reshape(real(Sq_est),Nx,Nz)';
        A=Sq_est;
		A(A<0)=0;

        figure(figSrc)
        subplot(1,2,2)
		imagesc(x_grid,z_grid,10*log10(A));
		%caxis([-30 5])
	    xlabel('x (m)')
	    ylabel('z (m)')
	    title(['MCMC FA, denoising = ' methods{i}])
	    hold on
	    plot(x_src,z_src,'*r')
	    colorbar
        
        borne=max(max(10*log10(A)));
	    caxis([borne-20 borne])
		set(gca,'Ydir','Normal')
%%
		    
		%%%-----------------------------------------------------------------------------
		%%% Save Data
		%%%-----------------------------------------------------------------------------
        
        %récupère les amplitudes à l'emplacement des vraies sources
        %{
 		for s=1:Nsrc
             amp_src(s) = BF_map(int32(ind_x_src(s)),int32(ind_z_src(s)));
        end
        
        
        %récupère les Nsrc plus grandes amplitudes de source
        amp_src = sort(BF_map(:),'descend');
        amp_src = amp_src(1:Nsrc)';
        %}
        
        %save Sc
        clear amp_src;
        amp_src(1,:,:)=Sq_est';
        
        h5write(file_name.map{i},'/SNR',SNR(j),[j 1],[1 1]);
        h5write(file_name.map{i},'/Amplitude_Source',amp_src,[j 1 1],[1 Nx Nz]);
        
        %Save BF map
        clear amp_src;
        amp_src(1,:,:)=BF_map;
        
        h5write(file_name.bf{i},'/SNR',SNR(j),[j 1],[1 1]);
        h5write(file_name.bf{i},'/Amplitude_Source',amp_src,[j 1 1],[1 Nx Nz]);
       	
        
        %Save Acoustic power
        %MCMC
        Gmap=GreenFreeField(r_map,r_map,freq,-fftsign,c);
        Gmap(logical(eye(Nx*Nz)))=0;
        Wmap(i,j)=1.2*c*k^2/(8*pi)*trace(real(Sq_est_corr)) + 1.2*c*k*0.5*real(trace(Gmap*1i*Sq_est_corr));
        %le resultat c-i-dessus est egal à 
        %r=1i*Gmap*4*pi/k;
        %r(logical(eye(Nmap)))=1;
        %r=r*(1.2*c*k*k/(8*pi));
        %trace(real(r)*Sq_est_corr)
        
        %et à Wmap(i,j)=1.2*c*k^2/(8*pi)*trace(real(Sq_est_corr)) - 1.2*c*k*0.5*imag(trace(Gmap*Sq_est_corr));
        

        %Reference power (real sources)
        Gsrc=GreenFreeField(r_src,r_src,freq,-fftsign,c);
        Gsrc(logical(eye(Nsrc)))=0;
        Wref(i,j)=1.2*c*k^2/(8*pi)*trace(real(Sq)) + 1.2*c*k*0.5*real(trace(Gsrc*1i*Sq));

        %beamforming
        Wbf(i,j) = 1.2*c*k^2/(8*pi)*trace(BF_map(int32(ind_x_src),int32(ind_z_src)));
        
        
        h5write(file_name.bf{i},'/Puissance',Wbf(i,j),[1 j],[1 1]);
        h5write(file_name.map{i},'/Puissance',Wmap(i,j),[1 j],[1 1]);
	end
end




%%%-----------------------------------------------------------------------------
%% % Display results
%%%-----------------------------------------------------------------------------
k=2*pi*freq/c;

figure
hold on
legend('show');
for i=1:length(methods)
   t=plot(SNR,real(10*log10(Wmap(i,:)./Wref(i,:))),'x--','DisplayName',['MCMC, ' methods{i}]);
   plot(SNR,real(10*log10(Wbf(i,:)./Wref(i,:))),'o-','DisplayName',['BF, ' methods{i}],'color',t.Color);

end
%plot(SNR,zeros(1,length(SNR)),':')
xlabel('SNR (dB)')
ylabel('$10 \log_{10}(\hat{W}/W_{ref})$')

%%
% Fidélité aux données
%{
Sy_est=G_map_mic*Sq_est_corr*G_map_mic' + mean(beta2_est(round(Nrun/2):end,:,1))*eye(Nmic);
figure
imagesc(10*log10(real(Sy./Sy_est)))

Sp_est=G_map_mic*Sq_est_corr*G_map_mic';
figure
imagesc(10*log10(real(Sp./Sp_est)))
%}

%%%-----------------------------------------------------------------------------
%% % 
%%%-----------------------------------------------------------------------------
%{

k=2*pi*freq/c;
%reconstructed power
Gmap=GreenFreeField(r_map,r_map,freq,-fftsign,c);
Gmap(logical(eye(Nx*Nz)))=0;
W=1.2*340*k^2/(8*pi)*trace(MAP.Sc) + 0.5* real(trace(Gmap*MAP.Sc));

%beamforming
Wbf = 1.2*340*k^2/(8*pi)*trace(BF_map(int32(ind_x_src),int32(ind_z_src)));

%real power
Gsrc=GreenFreeField(r_src,r_src,freq,-fftsign,c);
Gsrc(logical(eye(Nsrc)))=0;
Wref=1.2*340*k^2/(8*pi)*trace(real(Sq)) + 0.5* real(trace(Gsrc*Sq));
%}









