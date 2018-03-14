%Alice Dinsenmeyer; hiver 2017-2018
clear all
%close all

addpath('/home/adinsenmey/Bureau/Codes_Exterieurs/codes_debruitage_jerome')
addpath(genpath('/home/adinsenmey/Bureau/these_2017/Debruitage'))

freq = 3800;
Nsrc =1:5:95;
Mw=9000;
rho=0;
SNR = 10;
j=1;
for i=1:length(Nsrc)
i
    %%% Generate data
    [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw,SNR);
    %Sy(:,:,i)=Sp+diag(diag((Sn)));
    Rang2(i,j)=rank(Sp);

    CSM = Sy(:,:);
    d_ref(:,i,j)=real(diag(Sy-Sn));
    
    close all
    figure(55)
    plot(real(10*log10(sort(real(eig(Sn)),'descend'))))
    hold on
    plot(real(10*log10(sort(real(eig(Sp)),'descend'))));
    plot(real(10*log10(sort(real(eig(Sy)),'descend'))));
    plot(real(10*log10(sort(real(eig(Sy-Sn)),'descend'))));

    pause(1)


end


%save('cvx_rang','err_cvx','Rang2','d_cvx', 'd_ref');
%save('linprog_rang','err_linprog','Rang2','d_linprog','d_ref');
%save('AP_it_rang','err_it','Rang2','d_it','d_ref');




