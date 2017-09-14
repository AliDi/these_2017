function [Sqq_est , Sqq_diag , W ] = beamforming(Spp , G_map_mic)
%Estimates sources through beamfoming method. 
%
%Spp : (M x M ) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap): Transfert matrix from map points to mic position
%
%Returns :
%Sqq_est : (Nmap x Nmap ) matrix of source cross spectra
%Sqq_diag : (1 x Nmap ) diagonal terms of Sqq_est, ie source auto spectra
%W : (Nmap x M ) steering vectors 

	[M Nmap ]=size(G_map_mic);
        %%% Calculates weighting vectors (source scaling)
        L=diag(sum(abs(G_map_mic(:,:)).^2,1).^-1);%(G_map_mic(:,n,f)'*G_map_mic(:,n,f))^(-1); %scalar
            
        %%% Calculate sterring vectors
        W(:,:)=L*G_map_mic(:,:)';
        %W(n,:,f)=G_map_mic(:,n,f)'.*(G_map_mic(:,n,f)'*G_map_mic(:,n,f))^(-1);
		
		%%% Estimates sources  
		Sqq_est(:,:)=W(:,:)*Spp(:,:)*W(:,:)';
		Sqq_diag(:)=diag(Sqq_est(:,:)); %save auto spectra
end