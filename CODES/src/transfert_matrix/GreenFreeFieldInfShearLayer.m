function G=GreenFreeFieldInfShearLayer(r_src , r_mic , f , Gsign ,Mach , c , n_SL , pos_SL )
%G=GreenFreeFieldInfShearLayer(r_src , r_mic , f , Gsign ,Mach , c , n_SL , pos_SL )
%Calculates the acoustic transfert matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic
%through a infinitly thin shear layer (SL)
%n_SL is a vector normal to the interface flow/quiescent anechoic room
%pos_SL is the coordinate of a point of the planar interface
%c is the speed of sound in the room
%
%freq : excitation frequencies
%r_mic : (M x 3) matrix (column 1 : x, column 2 : y, column 3 : z)
%r_src : (N x 3) matrix (idem)
%G : (M x N) matrix
%
%/!!\ vector are signed !
%
    [M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    
    %reshape vectors
    Mach=reshape(Mach,1,3);
    n_SL=reshape(n_SL,1,3);
    pos_SL=reshape(pos_SL,1,3);
    
    r_mic=double(r_mic);
    r_src=double(r_src);
    
    Ra=pos_SL*n_SL'; %gives the distance from origine to shear layer
    
    k_out=2*pi*f/c; %wavenumber in the room
    
    for i=1:N
    
    	for m = 1:M
    		
    		u=r_mic(m,:)-r_src(i,:); %direction of propagation from source to receiver
    		u=u/norm(u); %make it unitary
    	
    		k_in = k_out/(1+Mach*u'); %wavenumber in the flow
    		
    		u_in = u+Mach; %direction of the propagation in the flow
    		u_in = u_in./norm(u_in); %make it unitary
    		
    		R_in = abs((Ra- r_src(i,:)*n_SL')/(u_in*n_SL')); %distance of propagation in the flow (from source to m)
    		
    		R_out = norm(r_mic(m,:)-r_src(i,:) - R_in*u_in);%distance of propagation in the room (from interface to mic)
    		

    		
    		G(m,i) = exp(Gsign*1j*k_in*R_in + Gsign*1j*k_out*R_out) /(4*pi*(R_in+R_out));
    	end
    end          
    		
end
