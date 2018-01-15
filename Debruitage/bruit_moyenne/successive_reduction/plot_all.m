load('AP_it_Mw.mat');
load('cvx_Mw.mat');
load('linprog_Mw.mat');


Nsrc =[20 60 80 96];
Mw=round(logspace(log10(10),log10(50000),50));
Mw2=round(logspace(log10(10),log10(50000),25));

figure

C = ['k','b','r'];
j=1;
for i=[ 1 3 4]

    semilogx(Mw,10*log10(err_cvx(i,:)),'color',C(j),'linestyle','--');
    hold on
    semilogx(Mw,10*log10(err_it(i,:)),'color',C(j),'linestyle',':');
    semilogx(Mw2,10*log10(err_linprog(i,:)),'color',C(j),'linestyle','-.');
   j=j+1; 
end   

%l=line([93 93],[-30 0],'color','black');
%l.Color(4)=0.3;
legend('Hald','Basic AP','Dougherty','location','southwest');
set(legend,'Box','off')
text(6*10^4, -22,'R=20')
text(6*10^4,-14,'R=80','color', C(2))
text(6*10^4,-11.5,'R=93','color', C(3))
xlim([Mw(1) 3*10^5])
ylim([-22.5 -7])
xlabel('Number of snapshots $M_w$')
ylabel('Relative error on diag(${S_p}$) (dB)')
plot_fig(gcf,11,8.5)
