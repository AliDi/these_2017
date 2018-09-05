close all
clear all
addpath('/home/adinsenmey/Bureau/these_2017/BeBeC2018_Debruitage/Generate_Spectra')
freq=15000;
Nsrc=20;
rho=0;
Mw=1221;
SNR=10;

[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR , 1221,SNR);
Sy(logical(eye(93)))=real(diag(Sy));
Sp(logical(eye(93)))=real(diag(Sp));

modif(:,:,1)= imag(Sy - Sp ) ./ imag(Sy) ;

%RPCA
Sp_den=h5read('../data/RPCA_Mw.h5','/Sp_real_Mw_1221') + 1i*h5read('../data/RPCA_Mw.h5','/Sp_imag_Mw_1221') ;

for i=1:101
    err(i)=norm(diag(Sp)-diag(Sp_den(:,:,i)),'fro')./norm(diag(Sp),'fro');
end
[err_min, ind_min]=min(err);
Sp_den=squeeze(Sp_den(:,:,ind_min));
Sp_den(logical(eye(93)))=real(diag(Sp_den));


modif(:,:,2)= imag(Sy - Sp_den ) ./ imag(Sy) ;


%MCMC
Sp_den=h5read('../data/MCMC_Mw.h5','/Sp_real_Mw_1221') + 1i*h5read('../data/MCMC_Mw.h5','/Sp_imag_Mw_1221') ;
Sp_den(logical(eye(93)))=real(diag(Sp_den));

modif(:,:,3)= imag(Sy - Sp_den ) ./ imag(Sy) ;

% Figure
seuil=5/100; %  5 pourcent par exemple
% for i=1:3
%     figure
%     spy(abs(modif(:,:,i))>seuil)
% end

for i=1:3
    figure
    imagesc(abs(modif(:,:,i)).*100)
end

clip=200;
cmap=viridis(64);

for i=1:3
    m=abs(modif(:,:,i))*100;
    m(m<=seuil*100)=0;


    m(m>clip)=clip;

    indx=find(m);
    [Nx Ny]=size(m);
    m=full(m(indx));
    ns = length(indx);
    [ix iy]=ind2sub([Nx Ny],indx);
    imap =m(:);


    figure, hold on
    colormap(cmap)
    scatter(iy,ix,[],imap,'Marker','.','SizeData',200)
    set(gca,'ydir','reverse')
    axis equal;
    xlabel(['nz = ' num2str(ns) '/8649'])
    axis([0 Nx 0 Ny])
    box on
    colorbar
    caxis([min(m) clip])
end





