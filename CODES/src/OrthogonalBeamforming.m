function [Sqq_diag , W , Nsrc] = OrthogonalBeamforming(Spp , G_map_mic , mic_weight)
%Estimates sources through orthogonal beamfoming method (Sarradj 2010)
%
%Spp : (M x M ) matrix of mic pressures cross spectra
%G_map_mic (M x Nmap): Transfert matrix from map points to mic position
%mic_weight : (1 x M) weight vector for microphones
%
%Returns :
%Sqq_diag : (Nmap x Nsrc) diagonal terms of Sqq_est, ie source auto spectra, for each source given by the highest eigenvalues of Spp
%W : (Nmap x M ) steering vectors 
%
       

	[M Nmap ]=size(G_map_mic);
	
	if (~exist('mic_weight'))
		mic_weight=ones(1,M);
	else
		mic_weight=reshape(mic_weight,1,M);
	end	
    
    [V Lambda] = eig(Spp);
    
    
    %%sort eigenvalues and eigenvectors 
	[Lambda,ind] = sort(diag(Lambda),'descend');
    V = V(:, ind);
    
    
    Lambda=diag(Lambda); 
    
    figure
    stem(diag(Lambda))
    title('Eigenvalues of $S_{pp}$')
    Nsrc = input('Number of sources ? ');
       
    Lambda = abs(real(Lambda));
    

    
    for i=1:Nsrc
    	Spp_i = V(:,i) * Lambda(i ,i ) * V(:,i)';
    	
    	for n=1:Nmap
            %%% Calculates weighting vectors (source scaling)
            L=sum(abs(G_map_mic(:,n)).^2,1).^-1; 
            
            %%% Calculate sterring vectors
            W(n,:)=L*G_map_mic(:,n)'.*mic_weight;
            
            %%% Estimates sources  
            Sqq_diag(n,i)=W(n,:)*Spp_i(:,:)*W(n,:)'; %save only auto spectra            
        end		
    	
    end 
end
