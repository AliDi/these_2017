function [cleanmap , dirtymap ] = clean_psf(LoopGain , Spp , G_map_mic , nb_itmax , trim , mic_weight)
%CLEAN-PSF (from Sijtsma, 2007)
%
%LoopGain (scalar 0< <1) : fraction of max value substracted too clean map
%Spp : (M x M ) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap ): Transfert matrix from map points to mic position
%nb_itmax (scalar) : max number of iterations
%trim : 'on' or 'off' to remove or not diagonal of Spp for each
%calculation
%mic_weight : (1xM) weight applied to each microphone
%
%cleanmap : (1xNmap) cleanmap
%dirtymap : (1xM) residuals

    [M Nmap]=size(G_map_mic);
    cleanmap=zeros(Nmap,1); 
    
	if (~exist('mic_weight'))
		mic_weight=ones(1,M);
	else
		mic_weight=reshape(mic_weight,1,M);
	end	


    %Initialize stop criterion
    count=1;
    StopCriterion=0;

    while ~StopCriterion

        count = count + 1

        %%% Trim Spp diagonal
        if strcmpi(trim,'on')
           Spp(logical(eye(M)))=0;               
        end 

        %%%  Obtain dirty map from beamforming
        [dirtymap , W] = beamforming(Spp(:,:) , G_map_mic(:,:) , mic_weight);

        %%% Find max of dirty and its position
        [valmax index_rmax] = max(real(dirtymap));
        valmax = valmax*LoopGain;

        %%% Find steering vector associated to rmax
        %wmax=W(index_rmax,:);

        %%% Update clean
        cleanmap(index_rmax) = cleanmap(index_rmax) + valmax;

        %%% Calculate Sppclean induced
        SppOld=Spp(:,:);
        Spp(:,:) = Spp(:,:) - valmax * G_map_mic(:,index_rmax)*G_map_mic(:,index_rmax)';

        %%%Stop criterion
        if ( norm(Spp(:,:),'fro') > norm(SppOld,'fro') || count > nb_itmax || valmax <0)
            StopCriterion=1;
            cleanmap(index_rmax) = cleanmap(index_rmax) - valmax;
            Spp(:,:) = Spp(:,:) + valmax * G_map_mic(:,index_rmax)*G_map_mic(:,index_rmax)';
        end

        %figure(5);
        %imagesc((real(cleanmap)));
        %colorbar
        %pause(1);
    end
    disp(['Number of iteration for CLEAN-PSF : ' num2str(count-1) ]);
   
end
