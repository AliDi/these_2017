clear all; close all;

load('s.mat');
rng(s);


jaune=[255 190 0]./255;
jauneclair=[255,223,128]./255;
rouge=[255 69 0]./255;
rougeclair=[255,162,128]./255;

M = 500;
meanX = 5.0; sigmaX = 10.0;
meanY = 0.0; sigmaY = 1;
rho=0.8;
x = zeros((M+1));
y = zeros((M+1));
x(1)=meanX;
y(1)=meanY;
c=sqrt(1-rho^2);

for i=2:(M+1)
    x(i)=meanX+sigmaX*rho*(y(i-1)-meanY)/sigmaY+ sigmaX*c*randn;
    y(i) = meanY+sigmaY*rho*(x(i)-meanX)/sigmaX+ sigmaY*c*randn;
end
delta = 0.025;
x1=-25:delta:35;
x2=-3:delta:3;
[X, Y] =meshgrid(x1,x2);
covariance=rho*sigmaX*sigmaY;

Z=mvnpdf([X(:) Y(:)],[meanX meanY],[sigmaX^2 covariance ; covariance sigmaY^2]);
Z=reshape(Z,size(X));


fig = figure(1);


z=subplot(2,2,2);
set(gca,'TickLabelInterpreter','latex')
set(gca,'DefaultTextInterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')
hold on
contour(X, Y, Z)
title('Gibbs Sampling')
grid on
xlim([min(x1) max(x1)])
ylim([min(x2) max(x2)])
box('on')
set(gca,'xticklabel',[],'yticklabel',[])

z2=subplot(2,2,1);
set(gca,'TickLabelInterpreter','latex')
set(gca,'DefaultTextInterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')
hold on
d=(max(x2)-min(x2))/round(sqrt(M));
plot(x2,normpdf(x2,meanY,sigmaY)*M*d,'k--')
h2=plot(0);
xlim([min(x2) max(x2)])
ylim([0 60])
p2=get(z2,'position');
set(z2,'position',[1.3*p2(1) p2(2)-0.4*p2(2) p2(3)-0.6*p2(3) p2(4)+0.5*p2(4)])
view([-90 90])
set(gca,'XTick',[-2  0 2])
xlabel('$\theta_1$')
box('on')
grid('on')



z3=subplot(2,2,4);
set(gca,'TickLabelInterpreter','latex')
set(gca,'DefaultTextInterpreter','latex')
set(gca,'DefaultLegendInterpreter','latex')
hold on
d=(max(x1)-min(x1))/round(sqrt(M));
plot(x1,normpdf(x1,meanX,sigmaX)*M*d,'k--')
h3=plot(0);
xlim([min(x1) max(x1)])
ylim([0 60])
p3=get(z3,'position');
set(z3,'position',[p3(1)-0.4*p3(1) 1.3*p3(2) p3(3)+0.5*p3(3) p3(4)/2])
xlabel('$\theta_2$')
box('on')
grid('on')



p=get(z,'position');
set(z,'position',[p(1)-0.4*p(1) p(2)-0.4*p(2) p(3)+0.5*p(3) p(4)+0.5*p(4) ])
pause(0.1)

print(fig,'-dpng','-r150',['img/' num2str(1,'%05.0f.png')])

%frame = getframe(fig);
%im = frame2im(frame);
%[imind,cm] = rgb2ind(im,256);
%filename='test.gif';
%imwrite(imind,cm,filename,'gif', 'DelayTime',0.005,'Loopcount',inf);


%%

index = 1;
for index=2:M
    
    subplot('position',z2.Position)
    delete(h2)
    h2=histfit(y(1:index),round(sqrt(M)));
    set(h2(1),'FaceColor',rougeclair,'EdgeColor',rouge);
    set(h2(2),'Color','k','LineWidth',1);
    
    
    subplot('position',z3.Position)
    delete(h3)
    h3=histfit(x(1:index),round(sqrt(M)));
    set(h3(1),'FaceColor',jauneclair,'EdgeColor',jaune);
    set(h3(2),'Color','k','LineWidth',1);
    set(gca,'XTick',[-20  0 20])
    
    subplot('position',z.Position)
    plot(x(index),y(index), '.k');
    xlim([min(x1) max(x1)])
    ylim([min(x2) max(x2)])
    
    %frame = getframe(fig);
    %im = frame2im(frame);
    %[imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    %imwrite(imind,cm,filename,'gif','DelayTime',0.005,'WriteMode','append');
    print(fig,'-dpng','-r150',['img/' num2str(index,'%05.0f.png')])

end

%positionne les valeurs opt
subplot('position',z2.Position)
set(gca,'XTick',[-2  meanY 2])
set(gca,'XTicklabel',{'-2'  '$\theta_1^*$' '2'})
line([meanY meanY],[0 60],'color',rouge,'linewidth',2)

subplot('position',z3.Position)
set(gca,'XTick',[-20  0 meanX 20])
set(gca,'XTicklabel',{'-20'  '0' '$\theta_2^*$' '20'})
line([meanX meanX],[0 60],'color',jaune,'linewidth',2)

subplot('position',z.Position)
line([min(x1) meanX],[meanY meanY],'color',rouge,'linewidth',2)
line([meanX meanX],[min(x2) meanY],'color',jaune,'linewidth',2)

print(fig,'-dpng','-r150',['img/' num2str(index+1,'%05.0f.png')])


%frame = getframe(fig);
%im = frame2im(frame);
%[imind,cm] = rgb2ind(im,256);

% Write to the GIF File
%imwrite(imind,cm,filename,'gif','DelayTime',2,'WriteMode','append');





