function G=GreenFreeField(r_src , r_mic , f , Gsign)
%Calculates the acoustic transfert matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic.
%
%freq : excitation frequencies
%r_mic : (M x 3) matrix (row 1 : x, row 2 : y, row 3 : z)
%r_src : (N x 3) matrix (idem)
%G : (M x N) matrix
%
    [M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    c=340; %speed of sound in air, m/sbeam
    k=2*pi*f/c;
    
    r_mic=double(r_mic);
    r_src=double(r_src);
    
    R=pdist2(r_mic , r_src);
    
    
    
    G = exp(Gsign*j*k*R)./(4*pi*R);          
end
