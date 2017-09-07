function [cleanmap , dirtymap ] = clean_psf(LoopGain , Spp , G_map_mic , nb_itmax)
%CLEAN-PSF (from Sijtsma, 2007)
%
%LoopGain (scalar 0< <1) : fraction of max value substracted too clean map
%Spp : (M x M x nb_F) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap x nb_F): Transfert matrix from map points to mic position
%F is supposed to be declared as global : (1 x nb_F) frequencies 
%nb_itmax (scalar) : max number of iterations
    global F;
    [M Nmap nb_F]=size(G_map_mic);
    cleanmap=zeros(Nmap,nb_F); 
    %Initialize stop criterium
    
    for iteration = 1:nb_itmax
        disp(['Iteration CLEAN : ' num2str(iteration) '/' num2str(nb_itmax)]);
        %%%  Obtain dirty map from beamforming
        [~ , dirtymap , W] = beamforming(Spp , G_map_mic);
        figure(5);
        imagesc((abs(dirtymap)));
        colorbar
        pause(1);
        for f=1:nb_F %CLEAN applied for each frequency indepently, for simplicity
            %%% If source map is 2D or 3D
            %dirty(:,:,f)=reshape(Sqq_diag(:,f),Nx,Ny,Nz);
            
            %%% Find max of dirty and its position
            [valmax index_rmax] = max(real(dirtymap(:,f)));
            index_rmax
            valmax
            valmax = valmax*LoopGain
            
            %%% Find steering vector associated to rmax
            wmax=W(index_rmax,:,f);
            
            %%% Update clean
            cleanmap(index_rmax,f) = cleanmap(index_rmax,f) + valmax;
            
            %%% Calculate Sppclean induced
            Spp(:,:,f) = Spp(:,:,f) - valmax * wmax' * wmax;
            
           %%% Trimm Spp
           tmp=Spp(:,:,f);
           tmp(logical(eye(M)))=0;
           Spp(:,:,f)=tmp;
           
        end
    end
    disp(['Number of iteration for CLEAN : ' num2str(iteration) ]);
end