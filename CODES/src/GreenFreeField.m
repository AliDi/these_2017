function G=GreenFreeField(r_src , r_mic)
%Calculates the acoustic transfert matrix for free field propagation (3D),
%from source points r_src to microphone points r_mic.
%
%freq : excitation frequencies
%r_mic : (M x 3) matrix (row 1 : x, row 2 : y, row 3 : z)
%r_src : (N x 3) matrix (idem)
%
%G : (M x N x nb_freq) matrix
    global F;
    [M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    c=343; %speed of sound in air, m/sbeam
    k=2*pi*F./c;
    for m=1:M
       for n=1:N
          R = norm(r_src(n,:)-r_mic(m,:));
          if R==0
              disp('/!\ Some source and mic points are superimposed.')
              return
          else
          G(m,n,:) = exp(-j*k*R)./(4*pi*R);          
          end
       end        
    end
end