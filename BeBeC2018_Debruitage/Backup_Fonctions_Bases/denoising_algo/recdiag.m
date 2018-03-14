function [matout,n,errec, diag_matout]=recdiag(mat,ro,nmax,tol,timetol)
%Code de Quentin Leclere pour le debruitage de la diagonale (2015)

if nargin==1
    nmax=30;
    tol=1e-6;
    timetol=30;
    ro=1;
end
if nargin==2
    tol=1e-6;
    timetol=30;
end

if nargin==3
    timetol=30;
end

matout = mat-diag(diag(mat));
continuer=1;
n=0;
tic
while continuer
    n=n+1;
    [V,D]=eig(matout);       
    d=diag(real(D));
    iok=find(d>=-tol);
    ipos=find(d>=0);
    matout = V(:,ipos)*diag(d(ipos))*V(:,ipos)';
    if length(iok)==length(d) | n>nmax | toc>timetol 
        continuer=0;
    else
        matout = (mat-diag(diag(mat)))+ro*diag(diag(matout));
    end
    diag_matout(:,n)=diag(matout);
end


errec=norm(matout-diag(diag(matout)) - mat+diag(diag(mat)), 'fro') ./ norm(mat-diag(diag(mat)), 'fro');


if nargout==0
    figure
    subplot(2,2,1);
    imagesc(abs(mat));
    clim = get(gca,'clim');
    subplot(2,2,2);
    imagesc(abs(matout));
    set(gca, 'clim', clim);
    title([num2str(round(errec*100)) '%'])
    
    subplot(2,2,3);
    plot(10*log10((real(diag(mat)))));
    hold on
    plot(10*log10((real(diag(matout)))), 'r');
    grid
    title(num2str(n))
    subplot(2,2,4);
    plot(sort(real(eig(mat)), 'descend'));
    hold on
    plot(sort(real(eig(matout)), 'descend'),'r');
    grid on
    
end

    
    
    
