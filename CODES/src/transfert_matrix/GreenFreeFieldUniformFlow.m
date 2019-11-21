function G=GreenFreeFieldUniformFlow(r_src , r_mic , f , Gsign , Mach , c )
%G=GreenFreeFieldUniformFlow(r_src , r_mic , f , Gsign , Mach c)
%Calculates the acoustic transfert matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic, SUBSONIC
%
%f : (scalar) excitation frequencie
%r_mic : (M x 3) matrix (column 1 : x, column 2 : y, column 3 : z)
%r_src : (N x 3) matrix (idem)
%G : (M x N) matrix
%Gsign : sign depending on convention for propagation 
%Mach : (1 x 3) Mach Vector
%
%From "NUMERICAL SIMULATION OF BEAMFORMING CORRECTION FOR DIPOLE SOURCE IDENTIFICATION", Yu Liu, 2008
%
    [M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    k=2*pi*f/c;
    
    Mach=reshape(Mach,1,3);
    
    %R=pdist2(r_mic , r_src); %distance between mic and src (1 x MN)
    
    beta2 = 1-norm(Mach)^2;
    
    %propagation vector betweeneach source and each receiver :
    R_vect = [ reshape(r_mic(:,1)-r_src(:,1)',M*N,1,1) , ...
    		   reshape(r_mic(:,2)-r_src(:,2)',M*N,1,1) , ...
    		   reshape(r_mic(:,3)-r_src(:,3)',M*N,1,1)]; 
    
    term_sqrt = sqrt( (sum(Mach.*R_vect,2)).^2 + beta2*sum(R_vect.^2,2) ); %intermediate calculation
    
    dr = (-sum(Mach.*R_vect,2) + term_sqrt )/beta2;%distance correction : dr/c = emmision time delay
  
        
    G = exp(Gsign*j*k*dr)./(4*pi*term_sqrt); 
	G=reshape(G,M,N);
         
end
