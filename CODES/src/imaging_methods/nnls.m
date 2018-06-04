function [Sqq] = nnls(Sqq_BF_diag , G_map_mic , W , options)
%[Sqq] = nnls(BF_Sqq_diag , G_map_mic , W )
%Deconvolution using Non negative least squares : 
%Sqq = min(|| Sqq_BF_diag - PSF*Sqq||²)
%
%--input--
%BF_Sqq_diag (1 x Nmap) : sources to be deconvolved, estimated with
%beamforming : source auto spectra
%G_map_mic (M x Nmap): Transfert matrix from map points to mic position
%W : (Nmap x M ) steering vectors
%Sqq : deconvolved source auto spectra
%

	[M Nmap ]=size(G_map_mic);   
    Sqq_BF_diag = real(reshape(Sqq_BF_diag, Nmap , 1));
    
    %%% Calculates PSF
    PSF = abs( W * G_map_mic ).^2;
    
    %%% NNLS
    [Sqq,~,~,exitflag,output]=lsqnonneg(PSF , Sqq_BF_diag,options);
    if ~exitflag
    	disp('!! nnls did not converg toward a solution !!')
    else
    	disp(['NNLS converged in ' num2str( output.iterations) ' iterations'])
    end
   
end
