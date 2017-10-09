function [Sqq]= damas(BF_Sqq_diag , G_map_mic , W , itermax)
%BF_Sqq_diag (1 x Nmap) : sources to be deconvolved, estimated with
%beamforming : source auto spectra
%G_map_mic (M x Nmap): Transfert matrix from map points to mic position
%W : (Nmap x M ) steering vectors
%Sqq : deconvolved source auto spectra

    [M Nmap ]=size(G_map_mic);   
    BF_Sqq_diag = real(BF_Sqq_diag);
    BF_Sqq_diag = reshape(BF_Sqq_diag,1 , Nmap);

    
    
    %%% Calculates PSF
    PSF = abs( W * G_map_mic ).^2;
    
    %%% Initialization
    Sqq = zeros(1 ,Nmap);
    Sqq_old = zeros(1 ,Nmap); %better convergence rate 
    %Sqq_old = BF_Sqq_diag;
    
    for i=1:itermax
    i
        for n=1:Nmap
            Sqq(n) = BF_Sqq_diag(n) - PSF(n , 1:n-1) * Sqq(1 : n-1 )' ...
                - PSF(n , n+1:end) * Sqq_old(n+1 : end)';
            
            %%% Ensure that power stays positive
            Sqq(n) = max(0,Sqq(n));
        end
        
        %%% Convergence criterion
        if max(abs(Sqq-Sqq_old))/mean(Sqq_old) < 0.05
            break;
        end
        
    Sqq_old = Sqq; 
    
    end
end
