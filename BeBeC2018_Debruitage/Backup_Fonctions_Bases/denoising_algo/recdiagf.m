function [matout]=recdiagf(mat)

invMfinez=pinvMfinez(size(mat,1));
matout=mat-diag(diag(mat)) + diag(exp(invMfinez*getbfinez(mat)));




if nargout==0
    figure
    subplot(2,2,1);
    imagesc(abs(mat));
    clim = get(gca,'clim');
    subplot(2,2,2);
    imagesc(abs(matout));
    set(gca, 'clim', clim);
    
    subplot(2,2,3);
    plot(10*log10((real(diag(mat)))));
    hold on
    plot(10*log10((real(diag(matout)))), 'r');
    grid
    subplot(2,2,4);
    plot(sort(real(eig(mat)), 'descend'));
    hold on
    plot(sort(real(eig(matout)), 'descend'),'r');
    grid on
    
end
