function [Sqq_diag , W ] = beamforming(Spp , G_map_mic , mic_weight)
%Estimates sources through beamfoming method. 
%
%Spp : (M x M ) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap): Transfert matrix from map points to mic position
%mic_weight : (1 x M) weight vector for microphones
%
%Returns :
%Sqq_diag : (1 x Nmap ) diagonal terms of Sqq_est, ie source auto spectra
%W : (Nmap x M ) steering vectors 
%
       

	[M Nmap ]=size(G_map_mic);
	
	if (~exist('mic_weight'))
		mic_weight=ones(1,M);
	else
		mic_weight=reshape(mic_weight,1,M);
	end	
    
    
        for n=1:Nmap
            %%% Calculates weighting vectors (source scaling)
            L=(sum(abs(G_map_mic(:,n)).^2,1)).^-1;%(G_map_mic(:,n,f)'*G_map_mic(:,n,f))^(-1); %scalar
            
            %%% Calculate sterring vectors
            W(n,:)=L*G_map_mic(:,n)'.*mic_weight;
            %W(n,:,f)=G_map_mic(:,n,f)'.*(G_map_mic(:,n,f)'*G_map_mic(:,n,f))^(-1);
            
            %%% Estimates sources  
            %Sqq_est(:,:)=W(:,:)*Spp(:,:)*W(:,:)';
            Sqq_diag(n)=W(n,:)*Spp(:,:)*W(n,:)'; %save only auto spectra            
        end		
end
