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
    
    
 
    
    figure
    stem(Lambda)
    title('Eigenvalues of $S_{pp}$')
    l(1)=Lambda(1);
    for i=2:M
		l(i)=l(i-1)+Lambda(i);
	end
	sum(Lambda)
    figure
    stem(l)
    ylabel('%')
    xlabel('N')    
    title('Énergie en % portée par la somme des N premières valeurs propres')
    
    Nsrc = input('Number of sources ? ');
       
    Lambda = abs(real(Lambda));
    Lambda=diag(Lambda);

    
    for i=1:Nsrc
    	Spp_i = V(:,i) * Lambda(i ,i ) * V(:,i)';
    	%figure
    	%imagesc(real(Spp_i))
    	
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
