function [matout,ii, x , diag_matout]=recdiagd(mat,nit,timetol)


matout = mat;
nm = size(mat,1);
[V,D]=eig(matout);
[d,ordre]=sort(diag(real(D)), 'descend');
V = V(:,ordre);
ii=0;
continuer=1;
sizeVold=1;
tic
while continuer
    ii=ii+1;
    A = abs(V').^2;
    for jj=sizeVold:size(V,2)
        b(jj)=abs(V(:,jj)'*mat*V(:,jj));
    end
    [t,x] = evalc('linprog(-1.*ones(nm,1),A,b,[],[],zeros(nm,1),real(diag(mat)));'); %evalc evite une sortie ecran de linprog
    %x=linprog(-1.*ones(nm,1),A,b,[],[],zeros(nm,1),real(diag(mat)));
    if size(x)==0
    	%matout=Inf*zeros(nm,nm);
    	break;
   	end
    matout = mat-diag(x);
    [Vn,Dn]=eig(matout);
    [dn,ordre]=sort(diag(real(Dn)), 'descend');
    Vn = Vn(:,ordre);
    sizeVold=size(V,2);
    V=[V Vn];
    
    elaps=toc;
    if  ii>nit | elaps>timetol
        continuer=0;
    end
    diag_matout(:,ii)=diag(matout);
end



%[V,D]=eig(matout);
%d =diag(D);
%iok = find(d>=0);
%matout = V(:,iok)*diag(d(iok))*V(:,iok)';
