function [q2psf, q2sc] = clean_sijtsma(R,k,S,nb,ro)



nm = size(R,1);
ns = size(R,2);

if isempty(k)
    G=R;
else
    G = exp(-1i*k*R)./(4*pi*R); % Green's functions
end


% w=zeros(ns,1);
%  for ii=1:ns
%      w(ii)=sum(sum(abs(G(:,ii)).^2 * abs(G(:,ii)').^2));
%  end
% W =diag(sqrt(w).^-1)* G' ;


%W=G';
    

%L = diag(sqrt(sum(abs(G).^2,1)).^-1); % source scaling terms 
%W = L*G'; % steering vectors


L = diag(sum(abs(G).^2,1).^-1); % source scaling terms 
Wq = L*G'; % steering vectors
%Wp=Wq;
Wp= diag(sqrt(diag(L)))*G';

%%%%%%%%%%%%%% CLEAN PSF %%%%%%%%%%%%%%
% initialize dirty CSM
D=S;
% initialize clean map
q2psf = zeros(ns,1); cc=0;
CONTINUE=1;
while CONTINUE
	% loop coutner
	cc=cc+1; 
	% search for the maximum in the dirty map 
    q2dirty=real(diag(Wp*D*Wp'));
    [~,imax]= findmax(q2dirty);  
    
    if ~isempty(imax)
        maxval=real(Wq(imax,:)*D*Wq(imax,:)');
        rocc=ro;
        %rocc=1-mean(q2dirty)./max(q2dirty);
        %rocc=max(rocc,ro);
        % add the maximum to the clean map (with a loop gain ro)
        q2psf(imax) = q2psf(imax) + rocc*maxval; 
        % save the cc-1 dirty map 
        Dold = D; 
        % remove the contribution of the maximum from the clean map
        D = D - rocc*maxval*G(:,imax)*G(:,imax)';
        % stoping criteria 
        %if norm(abs(D), 'fro')>norm(abs(D), 'fro') | cc>nb
        if norm(diag(D))>norm(diag(Dold)) | cc>nb  | maxval<0
            CONTINUE=0;
            q2psf(imax) = q2psf(imax) - rocc*maxval; 
        end
    else
        CONTINUE=0;
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%% CLEAN-SC %%%%%%%%%%%%%%
% initialize dirty CSM
D=S;
% initialize clean map
q2sc = zeros(ns,1); cc=0;
CONTINUE=1;
while CONTINUE
	% loop coutner
	cc=cc+1; 
	% search for the maximum in the dirty map 
    q2dirty=real(diag(Wp*D*Wp'));
    [~,imax]= findmax(q2dirty);
    
    if ~isempty(imax)
        maxval=real(Wq(imax,:)*D*Wq(imax,:)');
        rocc=ro;
        %rocc=1-mean(q2dirty)./max(q2dirty);
        %rocc=max(rocc,ro);
        % add the maximum to the clean map (with a loop gain ro)
        q2sc(imax) = q2sc(imax) + rocc*maxval;
        % calculate source component : D*W(imax,:)' are crossp between mics and max of the beam map
        h = D*Wq(imax,:)'./maxval;
        
        % save the cc-1 dirty map
        Dold = D;
        % remove (a part) of the coherent part of the microphone CSM
        D = D - rocc*maxval*h*h';
        % stoping criteria
        %if norm(abs(D), 'fro')>norm(abs(D), 'fro') | cc>nb
        if norm(diag(D))>norm(diag(Dold)) | cc>nb | maxval<0
            CONTINUE=0;
            q2sc(imax) = q2sc(imax) - rocc*maxval ;
        end
    else
        CONTINUE=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [maxval,imax]= findmax(input)

%[maxval,imax]= max(input);
imaxs = find(input(2:end-1)>input(1:end-2) & input(2:end-1)>input(3:end))+1;
[maxval,i] = max(input(imaxs));
imax = imaxs(i);




