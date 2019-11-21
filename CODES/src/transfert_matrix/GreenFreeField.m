function G=GreenFreeField(r_src , r_mic , f , Gsign ,c)
%G=GreenFreeField(r_src , r_mic , f , Gsign)
%Calculates the acoustic transfer matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic.
%
%freq : excitation frequency (scalar)
%r_mic : (M x 3) matrix (column 1 : x, column 2 : y, column 3 : z)
%r_src : (N x 3) matrix (idem)
%G : (M x N) matrix
%
    M= size(r_mic,1);
    N= size(r_src,1);
    k=2*pi*f/c;
    
    r_mic=double(r_mic);
    r_src=double(r_src);
    
    R=pdist2(r_mic , r_src);
    
    
    
    G = exp(Gsign*1i*k*R)./(4*pi*R);          
end
