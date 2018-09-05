addpath('/home/adinsenmey/Bureau/these_2017/BeBeC2018_Debruitage/Generate_Spectra')
data_path='./data/';

freq=15000;
rho=0;
SNR=10;
%=================================
%En fonction de Mw
%=================================
load([data_path 'AP_it_Mw.mat']); load([data_path 'cvx_Mw.mat']); load([data_path 'linprog_Mw.mat']);

Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),50));
Mw2=round(logspace(log10(10),log10(50000),25));

D{1}=d_cvx; D{2}=d_it; D{3}=d_linprog;

for j=1:length(Mw)
	for i=1:length(Nsrc)
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw(j),SNR);
		
        d=4;
        err(i,j,d)=norm(diag(Sp)-diag(Sy),'fro')/norm(diag(Sp),'fro');
    
        for d=1:2
		Sp_rec=Sy-diag(diag(Sy))+diag(D{d}(:,i,j));
		err(i,j,d)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
        end
	end
end	

d=3;
for j=1:length(Mw2)
	for i=1:length(Nsrc)
        
		%%% Generate data
		[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw2(j),SNR);

		Sp_rec=Sy-diag(diag(Sy))+diag(D{d}(:,i,j));
		err(i,j,d)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');	
	end
end

%%
figure
C = ['k','b','r'];
L={ {'--'} ,  {':'} , {'-.'} ,  {'-' , 'color' , [0,0,0,0.5]}};
j=1;
for i=[ 1 3 4]
	for d=[4 1 2]
        semilogx(Mw,10*log10(err(i,:,d)),'color',C(j),'linestyle',L{d}{:});
        hold on
    end
    d=3;
    semilogx(Mw2,10*log10(err(i,1:25,d)),'color',C(j),'linestyle',L{d}{:});
    
	j=j+1; 
end
%%	

%legend('Hald','Basic AP','Dougherty','location','southwest');
%set(legend,'Box','off')
text(20, -20,'K=20'); text(20,-18,'K=80','color', C(2)); text(20,-16,'K=93','color', C(3))
xlim([Mw(1) Mw(end)]); ylim([-21 -3]);
xlabel('Number of snapshots $M_w$'); ylabel('Relative error on ${S_p}$ (dB)');
plot_fig(gcf,11,8.5)


%%
clear err;
Mw=10^4;
SNR=10;
%=================================
%En fonction de Nsrc
%=================================
load([data_path 'AP_it_rang.mat']); load([data_path 'linprog_rang.mat']);
Nsrc2=Nsrc;
load([data_path 'cvx_rang.mat']);

D{1}=d_cvx; D{2}=d_it; D{3}=d_linprog;

j=1;
for i=1:length(Nsrc)
    %%% Generate data
    [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw,SNR);
    
    d=4;
    err(i,d)=norm(diag(Sy)-diag(Sp),'fro')/norm(diag(Sp),'fro');
    for d=1:3
        Sp_rec=Sy-diag(diag(Sy))+diag(D{d}(:,i,j));
        err(i,d)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
    end
end	

%%
figure
%hold on
C = ['k','b','r'];
L={ {'--'} {':'} {'-.'} {'-' , 'color' , [0,0,0,0.5]}};
j=1;

for d=[1 2 4]
    plot(Nsrc,10*log10(err(:,d)),'Color','k','linestyle',L{d}{:});
    hold on
    
end
d=3;
plot(Nsrc2,10*log10(err(1:length(Nsrc2),d)),'linestyle',L{d}{:},'Color','k');

%%	
%legend('Hald','Basic AP','Dougherty','location','southwest');
%set(legend,'Box','off')

xlim([Nsrc(1) Nsrc(end)])
ylim([-22 -12])
xlabel('Rank of $S_{pp}$')
ylabel('Relative error on ${S_p}$ (dB)')
plot_fig(gcf,11,8.5)
%%


%%
clear err;
Mw=10^4;
%=================================
%En fonction de SNR
%=================================
load([data_path 'AP_it_SNR.mat']); load([data_path 'cvx_SNR.mat']); load([data_path 'linprog_SNR.mat']);
Nsrc=20;

D{1}=d_cvx; D{2}=d_it; D{3}=d_linprog;


for i=1:length(SNR)
    %%% Generate data
    [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw,SNR(i));
    
    d=4;
    err(i,d)=norm(diag(Sp)-diag(Sy),'fro')/norm(Sp,'fro');

    for d=1:3
    Sp_rec=Sy-diag(diag(Sy))+diag(D{d}(:,i,j));
    err(i,d)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
    end
	
end	

%%
figure
hold on
C = ['k','b','r'];
L={ {'--'} {':'} {'-.'} {'-' , 'color' , [0,0,0,0.5]}};
j=1;

for d=1:4
    plot(SNR,10*log10(err(:,d)),'Color','k','linestyle',L{d}{:});    
end
%%	

xlim([SNR(1) SNR(end)])
ylim([-18 6])
xlabel('SNR (dB)')
ylabel('Relative error on ${S_p}$ (dB)')
plot_fig(gcf,11,8.5)
%%





