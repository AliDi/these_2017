% trace la pression des sources repropagées sur une sphère de de rayon ray.
figure
plot3(r_mic(:,1),r_mic(:,2),r_mic(:,3),'ob')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
plot3(r_src(:,1),r_src(:,2),r_src(:,3),'r*')
grid on

N=70;
theta=linspace(0,2*pi,N);
phi=linspace(-pi/2,pi/2,N);
ray=1;
[theta phi]=meshgrid(theta, phi);
x=ray.*cos(theta).*cos(phi);
y=ray.*sin(theta).*cos(phi);
z=ray.*sin(phi);


%figure
%plot3(x(:),y(:),z(:),'ob')
r_sp=[x(:) y(:) z(:)];

%directivité des sources reconstruites
Gsp=GreenFreeField(r_map,r_sp,freq,-fftsign,c);
Spp=(1.2*c*k)^2*Gsp*Sq_est_corr*Gsp';

ray=real(reshape(diag(Spp),N,N));

x=ray.*cos(theta).*cos(phi);
y=ray.*sin(theta).*cos(phi);
z=ray.*sin(phi);

h=figure;
surf(x,y,z,ray);
xlabel('x')
ylabel('y')
zlabel('z')
title('Reconstructed')





%reference
Gsp=GreenFreeField(r_src,r_sp,freq,-fftsign,c);
Spp_ref=(1.2*c*k)^2*Gsp*Sq*Gsp';

ray_ref=real(reshape(diag(Spp_ref),N,N));

x=ray_ref.*cos(theta).*cos(phi);
y=ray_ref.*sin(theta).*cos(phi);
z=ray_ref.*sin(phi);

figure;
surf(x,y,z,ray_ref);
xlabel('x')
ylabel('y')
zlabel('z')
xlim([min(min([x y z])) max(max([x y z]))])
ylim([min(min([x y z])) max(max([x y z]))])
zlim([min(min([x y z])) max(max([x y z]))])
title('Reference')
caxis([min(min(ray_ref)) max(max(ray_ref))])
colorbar;

figure(h)
xlim([min(min([x y z])) max(max([x y z]))])
ylim([min(min([x y z])) max(max([x y z]))])
zlim([min(min([x y z])) max(max([x y z]))])
caxis([min(min(ray)) max(max(ray))])
colorbar



