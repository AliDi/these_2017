clear all; close all;

M = 500;
meanX = 0.0; sigmaX = 10.0;
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
x1=-30:delta:30;
x2=-3:delta:3;
[X, Y] =meshgrid(x1,x2);
covariance=rho*sigmaX*sigmaY;

Z=mvnpdf([X(:) Y(:)],[meanX meanY],[sigmaX^2 covariance ; covariance sigmaY^2]);
Z=reshape(Z,size(X));


fig = figure(1);
z=subplot(2,2,2);
hold on
contour(X, Y, Z)
title('Gibbs Sampling')
grid on
xlim([-30 30])
ylim([-3 3])

z2=subplot(2,2,1);
plot(0)
xlim([-3 3])
ylim([0 60])
p2=get(z2,'position');
set(z2,'position',[p2(1) p2(2)-0.4*p2(2) p2(3)-0.6*p2(3) p2(4)+0.5*p2(4)])

z3=subplot(2,2,4);
plot(0)
xlim([-30 30])
ylim([0 60])
p3=get(z3,'position');
set(z3,'position',[p3(1)-0.4*p3(1) p3(2) p3(3)+0.5*p3(3) p3(4)/2])


p=get(z,'position');
set(z,'position',[p(1)-0.4*p(1) p(2)-0.4*p(2) p(3)+0.5*p(3) p(4)+0.5*p(4) ])

frame = getframe(fig);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
filename='test.gif';
imwrite(imind,cm,filename,'gif', 'DelayTime',0.05,'Loopcount',inf);

%%

index = 1;
for index=2:M

    subplot('position',z2.Position)
    h2=histfit(y(1:index),12);
    xlim([-3 3])
    ylim([0 60])
    set(h2(1),'FaceColor','green','EdgeColor','green');
    set(h2(2),'Color','k','LineWidth',1);
    view([-90 90])
    set(gca,'XTick',[-2  0 2])


    subplot('position',z3.Position)
    h2=histfit(x(1:index));
    xlim([-30 30])
    ylim([0 60])
    set(h2(1),'FaceColor','red','EdgeColor','red');
    set(h2(2),'Color','k','LineWidth',1);
    set(gca,'XTick',[-20  0 20])

    subplot('position',z.Position)
    plot(x(index),y(index), 'ob');
    xlim([-30 30])
    ylim([-3 3])
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File    
    imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append');
    
    
    
    
end