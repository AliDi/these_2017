function [cleanmap , dirtymap , count ] = clean_sc(LoopGain , Spp , G_map_mic , nb_itmax , trim)
%CLEAN-SC (from Sijtsma, 2007)
%
%LoopGain (scalar 0< <1) : fraction of max value substracted too clean map
%Spp : (M x M x nb_F) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap x nb_F): Transfert matrix from map points to mic position
%nb_itmax (scalar) : max number of iterations
%trim : 'on' or 'off' to remove or not diagonal of Spp for each
%calculation

    [M Nmap]=size(G_map_mic);
    cleanmap=zeros(Nmap,1); 
    
    %Initialize stop criterion
    count=1;
    StopCriterion=0;

    while ~StopCriterion

        count = count + 1;

        %%% Trim Spp diagonal
        if strcmpi(trim,'on')
           Spp(logical(eye(M)))=0;               
        end 

        %%%  Obtain dirty map from beamforming
        [dirtymap , W] = beamforming(Spp(:,:) , G_map_mic(:,:));

        %%% If source map is 2D or 3D
        %dirty(:,:,f)=reshape(Sqq_diag(:,f),Nx,Ny,Nz);

        %%% Find max of dirty and its position
        [valmax index_rmax] = max(real(dirtymap));
        valmax = valmax*LoopGain;

        %%% Find steering vector associated to rmax
        wmax=W(index_rmax,:);
        wmax=wmax'; %make wmax a colonn vector

        %%% Update clean
        cleanmap(index_rmax) = cleanmap(index_rmax) + valmax;
        
        %%% Calculate coherent sources propagator
        if strcmpi(trim,'on')
            h = G_map_mic(:,index_rmax);
            for iterH = 1:1000
                hOldValue = h;
                H = h*h';        
                H(~logical(eye(M))) = 0;            
                h = 1/sqrt(1+wmax'*H*wmax)*(Spp*wmax/valmax + H*wmax);
                if norm(h-hOldValue) < 1e-10
                    break;
                end
            end
        else
             h=Spp*wmax./valmax;
        end
        
        %%% Calculate Sppclean induced
        SppOld=Spp(:,:);
        Spp(:,:) = Spp(:,:) - valmax * h * h';

        %%%Stop criterion
        if ( norm(Spp(:,:),'fro') > norm(SppOld,'fro') || count > nb_itmax || valmax <0)
            StopCriterion=1;
            cleanmap(index_rmax) = cleanmap(index_rmax) - valmax;
            Spp(:,:) = Spp(:,:) + valmax * h *h';
        end

        %figure(5);
        %imagesc((real(cleanmap)));
        %colorbar
        %pause(1);
    end
    count=count-1;
    disp(['Number of iteration for CLEAN-SC : ' num2str(count) ]);
end
