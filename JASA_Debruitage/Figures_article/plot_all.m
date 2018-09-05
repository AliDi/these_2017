%close all;
clear all;

addpath('/home/adinsenmey/Bureau/these_2017/JASA_Debruitage/Data/Generate_Spectra')

lambda=0:0.01:1;
freq=15000;
rho=0;
%%
%param(1).method='cvx';
%param(1).name=[{'Mw'} {'Nsrc'} {'SNR'}];
%param(1).value=[ {round(logspace(log10(10),log10(50000),50)) } ... %Mw
%    {1:93} ... 											%Nsrc
%    {-10:10}  ];					%SNR
%    
%param(2).method='RPCA';
%param(2).name=[{'Mw'} {'Nsrc'} {'SNR'}];
%param(2).value=[ {round(logspace(log10(10),log10(50000),40)) } ... %Mw
%    {1:93} ... 											%Nsrc
%    {-10:10}  ]; 										%SNR
    
			
    
    
% Export erreur (en .mat) pour RPCA, MCMC et la référence non débruitée
%%====================================================

%Mw
%%
%{
Nsrc=[20 60 80 96];
SNR=10;
Mw=round(logspace(log10(10),log10(50000),40));
for j=1:length(Nsrc)
    for i=1:length(Mw)
        [Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(j) , rho , SNR , Mw(i),SNR);

        %NONE (référence sans débruitage)
        err_none(j,i)=norm(diag(Sp)-diag(Sy),'fro')/norm(diag(Sp),'fro');

%         %RPCA
%         Sp_rec=h5read('./data/RPCA_Mw.h5',['/Sp_real_Mw_' num2str(Mw(i))]) + 1i*h5read('./data/RPCA_Mw.h5',['/Sp_imag_Mw_' num2str(Mw(i))]);
%         for j=1:length(lambda)
%             err(i,j)=norm(diag(Sp)-diag(Sp_rec(:,:,j)),'fro')/norm(diag(Sp),'fro');
%         end
%         err_01=err(:,11);
%         err_opt=min(err,[],2);
%         err_RPCA= [err_01 err_opt];
% 
%         %MCMC
%         Sp_rec=h5read('./data/MCMC_Mw.h5',['/Sp_real_Mw_' num2str(Mw(i))]) + 1i*h5read('./data/MCMC_Mw.h5',['/Sp_imag_Mw_' num2str(Mw(i))]);
%         err_MCMC(i)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
    end
end
save('none_Mw','err_none','Mw','Nsrc');
%save('RPCA_Mw','err_RPCA','Mw');
%save('MCMC_Mw','err_MCMC','Mw');
%%

%Nsrc
Nsrc=1:93;
SNR=10;
Mw=10^4;
for i=1:length(Nsrc)
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc(i) , rho , SNR , Mw,SNR);
    
    %NONE (référence sans débruitage)
    err_none(i)=norm(diag(Sp)-diag(Sy),'fro')/norm(diag(Sp),'fro');
	
	
	%RPCA
	Sp_rec=h5read('./data/RPCA_Nsrc.h5',['/Sp_real_Nsrc_' num2str(Nsrc(i))]) + 1i*h5read('./data/RPCA_Nsrc.h5',['/Sp_imag_Nsrc_' num2str(Nsrc(i))]);
	for j=1:length(lambda)
		err(i,j)=norm(diag(Sp)-diag(Sp_rec(:,:,j)),'fro')/norm(diag(Sp),'fro');
	end
	err_01=err(:,11);
	err_opt=min(err,[],2);
	err_RPCA= [err_01 err_opt];
    
    %MCMC
    Sp_rec=h5read('./data/MCMC_Nsrc.h5',['/Sp_real_Nsrc_' num2str(Nsrc(i))]) + 1i*h5read('./data/MCMC_Nsrc.h5',['/Sp_imag_Nsrc_' num2str(Nsrc(i))]);
    err_MCMC(i)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
end
%save('none_Nsrc','err_none','Nsrc');
%save('RPCA_Nsrc','err_RPCA','Nsrc');
save('MCMC_Nsrc','err_MCMC','Nsrc');


%%


%SNR
Nsrc=20;
Mw=10^4;
SNR=-10:10;
for i=1:length(SNR)
	[Sq Sy Sp Sn] = generate_Spp_signal(freq, Nsrc , rho , SNR(i) , Mw,SNR(i));
    
    %NONE (référence sans débruitage)
    err_none(i)=norm(diag(Sp)-diag(Sy),'fro')/norm(diag(Sp),'fro');
	
	
	%RPCA
	Sp_rec=h5read('./data/RPCA_SNR.h5',['/Sp_real_SNR_' num2str(SNR(i))]) + 1i*h5read('./data/RPCA_SNR.h5',['/Sp_imag_SNR_' num2str(SNR(i))]);
	for j=1:length(lambda)
		err(i,j)=norm(diag(Sp)-diag(Sp_rec(:,:,j)),'fro')/norm(diag(Sp),'fro');
	end
	err_01=err(:,11);
	err_opt=min(err,[],2);
	err_RPCA= [err_01 err_opt];
    
    %MCMC
    Sp_rec=h5read('./data/MCMC_SNR.h5',['/Sp_real_SNR_' num2str(SNR(i))]) + 1i*h5read('./data/MCMC_SNR.h5',['/Sp_imag_SNR_' num2str(SNR(i))]);
    err_MCMC(i)=norm(diag(Sp)-diag(Sp_rec),'fro')/norm(diag(Sp),'fro');
end
%save('none_SNR','err_none','SNR');
%save('RPCA_SNR','err_RPCA','SNR');
save('MCMC_SNR','err_MCMC','SNR');

%}
%%
%% Diagonale reconstrution
close all
path='./data/';
param=[{'Mw'} {'Nsrc'} {'SNR'}];
method=[  {'cvx'} {'AP_it'} {'linprog'} {'none'}];
L={  {'-'} {'--'}   , {':'}, {'-', 'color' , [0.5 0.5 0.5]} };
M={{'none'} {'none'} {'none'} {'none'}};   %{ {'x'} ,  {'*'} , {'s'} ,  {'none'}};

for p=1:3
    for m=1:4
	

		load([path method{m} '_' param{p}])
        
		figure(p)
        hold on
        if strcmp(param{p},'Mw')
            set(gca, 'ColorOrder', [0    0.4470    0.7410 ; 0.8500    0.3250    0.0980 ; 0.9290    0.6940    0.1250 ], 'NextPlot', 'replacechildren');
            hold on
            plot(eval(param{p}),10*log10(eval(['err_' method{m} '([1 3 4],:)'])),'linestyle',L{m}{:},'Marker',M{m}{:}); %
        else
            set(gca, 'ColorOrder', [0    0 0  ], 'NextPlot', 'replacechildren');
            hold on
            plot(eval(param{p}),10*log10(eval(['err_' method{m} ])),'linestyle',L{m}{:},'Marker',M{m}{:}); %
        end
        xlim([min(eval(param{p})) max(eval(param{p}))])

		
    end
%     m=4;
%     figure(p)
%     load([path method{m} '_' param{p}])
%     plot(eval(param{p}),10*log10(eval(['err_' method{m}])),'linestyle',L{m}{:},'Marker',M{m}{:});
     legend('Convex','AP','Lin Prog','none')
end

figure(1)
set(gca,'xscale','log')
line([93 93],[-30 0],'linewidth',1,'color',[0.5 0.5 0.5])
legend('K=20','K=80','K=93')



%%

%{    
%%% Comparaison toutes les methodes
%%%====================================================
path='./data/';
param=[{'Mw'} {'Nsrc'} {'SNR'}];
method=[{'none'} {'cvx'} {'RPCA'}  {'MCMC'} {'EM'}];
%M={ {'none'} ,  {'x'} , {'s'} ,   {'d'}}; %{'*'}
L={ {'-', 'color' , [0.5 0.5 0.5]} {'-', 'color' ,'k'} {'--'} ,  {'--', 'color' ,'k'}  , {':', 'color' ,'k'} };


for p=1:3
    for m=1:length(method)
		load([path method{m} '_' param{p}])
        
		figure(p+10)
        set(gca, 'ColorOrder', [0 0 0 ;0.5 0.5 0.5], 'NextPlot', 'replacechildren');
		hold on
        if strcmp(method{m},'cvx')
      		plot(eval(param{p}),10*log10(eval(['err_' method{m} '(1,:)'])),'linestyle',L{m}{:},'linewidth',2);
        else
            plot(eval(param{p}),10*log10(eval(['err_' method{m}])),'linestyle',L{m}{:},'linewidth',2);
        end
		
    end
    xlim([min(eval(param{p})) max(eval(param{p}))])
    legend('Reference','DRec','RPCA, l=0.1','RPCA, l_{opt}','MCMC','EM')
    hold off
end

%finitions
figure(1+10)
set(gca,'xscale','log')
line([93 93],[-30 0],'linewidth',1,'color',[0.5 0.5 0.5])
legend('Reference','DRec','RPCA, l=0.1','RPCA, l_{opt}','MCMC','EM')


%% Choix de K pour MCMC/EM

method=[ {'MCMC'} {'EM'}];

for m=1:2
    
    load([path method{m} '_K'])
    figure
    hold on
    plot(K,10*log10(eval(['err_' method{m}])),'o-')
end


%}






%    
%%% Trace figure
%%%====================================================

%for p =1:2
%    for n=1:length(param(p).name) %pour chaque paramètre
%        
%        clear err
%        for v=1:length(param(p).value{n}) %pour chaque valeur du paramètre
%            err(v,:) = 	h5read(['./data/' param(p).method '_' param(p).name{n} '.h5'],['/' param(p).name{n} '_' num2str(param(p).value{n}(v))]);
%        end
%        
%        if param(p).method=='RPCA'
%            err_opt=min(err,[],2);
%            err_01=err(:,11);
%            
%            figure(n+10)
%            plot(param(p).value{n},10*log10(err_01))
%            hold on
%            plot(param(p).value{n},10*log10(err_opt))
%        else
%            figure(n+10)
%            hold on
%            plot(param(p).value{n},10*log10(err))
%        end
%    end
%end


%%% Réglage affichage
%figure(1+10)
%set(gca,'xscale','log')





