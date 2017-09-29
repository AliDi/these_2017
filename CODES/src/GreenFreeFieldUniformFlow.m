function G=GreenFreeFieldUniformFlow(r_src , r_mic , f , Gsign , Mach)
%Calculates the acoustic transfert matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic.
%
%freq : excitation frequencies
%r_mic : (M x 3) matrix (row 1 : x, row 2 : y, row 3 : z)
%r_src : (N x 3) matrix (idem)
%G : (M x N) matrix
%Gsign : sign depending on convention for propagation 
%Mach : (scalar) Mach Number
%FlowDirection  : scalar = 1,2 or 3 : number corresponding to the direction of the uniform flow.
%
%
%
    [M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    c=340; %speed of sound in air, m/sbeam
    k=2*pi*f/c;
    
    R=pdist2(r_mic , r_src);

	Rtilde = sqrt( pdist2(r_mic(:,1),r_src(:,1)).^2 +  )

	    
    
    G = exp(Gsign*j*k*R)./(4*pi*R); 
         
end