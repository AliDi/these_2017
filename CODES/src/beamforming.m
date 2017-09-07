function [Sqq_est , Sqq_diag , W ] = beamforming(Spp , G_map_mic)
%Estimates sources through beamfoming method. 
%
%Spp : (M x M x nb_F) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap x nb_F): Transfert matrix from map points to mic position
%F is supposed to be declared as global : (1 x nb_F) frequencies 
%
%Returns :
%Sqq_est : (Nmap x Nmap x nb_F) matrix of source cross spectra
%Sqq_diag : (Nmap x nb_F) diagonal terms of Sqq_est, ie source auto spectra
%W : (Nmap x M x nb_F) steering vectors 

	[M Nmap nb_F]=size(G_map_mic);
	for f=1:nb_F
	
		%%% Calculate weights/source scaling
		for m=1:M
			Lvect(m)=1/norm(G_map_mic(m,:,f));
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
end