function [Gden flag] = CanonicalCoherences(G,option)
%
% Canonical-coherence-based denoising
%
% Gden = CANONICALCOHERENCES(G,option)
%	G : CSM matrix (MxM) to denoise
%	Gden : denoised CSM
%
%	If option.threshold = 1, the thresholding of the canonical components is applied, as given by equation (18a) and (18b) in [1]
%
%	If option.iterative = 1, the iterative procedure is applied as described in [1]
%
% [1] Denoising of cross-spectral matrices using canonical coherence, J. Hald, JASA 2019
%
% Implemented by Alice Dinsenmeyer, 26 Sept. 2019
%
M=size(G,1);
if mod(M,2) || mod(M,6)
 %   error('The microphone number must divisble by 6');
end

% Selection of the subgroups
if option.iterative ~= 1
	
	indsub1(:,1) = 1:2:M; %indices for the first subgroup
	indsub2(:,1) = 2:2:M; %indices for the second subgroup
	
else
	
	%Preallocation
	indsub1=zeros(floor(M/2),3);
	indsub2=zeros(ceil(M/2),3);
	
	%first iteration
	indsub1(:,1) = 2:2:M; %indices for the first subgroup
	indsub2(:,1) = 1:2:M; %indices for the second subgroup
	
	%second iteration
	a=1:12:M;
	b=[1 3 5 8 10 12];
	ind= vec(a+b.'-1);
    ind=ind(1:floor(M/2));
	indsub1(:,2) = ind(ind<=M); %indices for the first subgroup
	indsub2(:,2) = setxor(1:M,indsub1(:,2)); %indices for the second subgroup
	
	%third iteration
	indsub1(:,3) = sort(randperm(M,floor(M/2))); %indices for the first subgroup
	indsub2(:,3) = setxor(1:M,indsub1(:,3)); %indices for the second subgroup

end

Gin = G;
Gout_old = 0;
for i = 1:size(indsub1,2)

	Gxx = Gin(indsub1(:,i),indsub1(:,i));
	Gyy = Gin(indsub2(:,i),indsub2(:,i));
    
    Gxy = Gin(indsub1(:,i),indsub2(:,i));
    [I J]=size(Gxy);
	
	Gxx_sqrt = sqrtm(Gxx);
	Gyy_sqrt = sqrtm(Gyy);
	
	K = inv(Gxx_sqrt) * Gxy * inv(Gyy_sqrt); %from eq (5)
	
	[U,Sigma,V] = svd(K);
    
    if option.threshold == 1 %thresholding
       if i==1
           sig_thresh = max(0.3*max(diag(Sigma)),0.2); %eq (18a)
       else
           sig_thresh = max(0.6*max(diag(Sigma)),0.2); %eq (18b). these thresholds works only for medium SNR
       end
       Sigma(Sigma < sig_thresh) = 0;           
    end   
    
	Sigma_sqrt = sqrt(Sigma); 
	
	P = Gxx_sqrt* U * Sigma_sqrt; %from eq (12a)
	Q = Gyy_sqrt * V * Sigma_sqrt';%from eq. (12b)
    
    P=P(:,1:min(I,J)); 
    Q=Q(:,1:min(I,J)); 
    PQ = P*Q';
    
    Gout = [P*P' PQ ; PQ' Q*Q'];
    Gout = unshuffle(Gout,[indsub1(:,i) ; indsub2(:,i)],[indsub1(:,i) ;indsub2(:,i)]); %reshuffle the rows of Gout
    
    Gout = Gout_old + Gout; %cumulated output, eq (17a)
    Gin = G-Gout; %next input
    Gout_old = Gout;
    
    %empirical factors from eq. (19a,b)
    if i==1
        alpha=trace(Gout)/trace(G);
    end
    if alpha < 0.25
        beta=1;
    else
        if i==1
            beta=0.19+0.133*(1-alpha);            
        elseif i==2
            beta=0.35+0.133*(1-alpha);
        end
    end
    
    s=eig(Gout);
    s=10*log10(s./max(s));
    Lambda = sum(s>-10); %number of eigenvalues in Gout that are no more than 10 db below the highest.
    
    if (i==1 & Lambda<= beta*M)
        break;
    elseif (i==2 & Lambda<= beta*M)
        break;
    end 
    
    
end
flag.count=i;
Gden = Gout;




