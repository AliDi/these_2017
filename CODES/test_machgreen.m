clear all;
close all;

addpath('~/Bureau/these_2017/CODES/src')

%%% Microphones position
array_length = 0.7; %in meters
M = 32; %number of microphones
z_mic = 0 * ones(1,M);

x_mic = linspace(0,0.7,M);
y_mic = zeros(1,M);
r_mic = [x_mic' y_mic' z_mic'];

%%% Define grid for source map
Nx=10;
Ny=10;
Nz=1;
Nmap = Nx * Ny * Nz; %number of source point on the constructed map

x_grid = linspace(-0.5,0.5,Nx);
y_grid = linspace(-0.5,0.5,Ny);
z_grid = linspace(0.7 , 0.7 , Nz);

[Xmap,Ymap, Zmap] = meshgrid(x_grid,y_grid,z_grid);
Xmap=Xmap';
Ymap=Ymap';
Zmap=Zmap';

r_src_map = [Xmap(:) Ymap(:) Zmap(:)]; % vector of estimate source point positions (map)
r_src_map = double(r_src_map);

r_src=double(r_src_map);
f = 100;
Gsign=1;
Mach=[0.1 0 0];


	%function
	[M, ~]= size(r_mic);
    [N, ~]= size(r_src);
    c=340; %speed of sound in air, m/sbeam
    k=2*pi*f/c;
    %R=pdist2(r_mic , r_src); %distance between mic and src (1 x MN)
    
    beta2 = 1-norm(Mach)^2;
    
    %propagation vector betweeneach source and each receiver :
    R_vect = [ reshape(r_mic(:,1)-r_src(:,1)',M*N,1,1) , ...
    		   reshape(r_mic(:,2)-r_src(:,2)',M*N,1,1) , ...
    		   reshape(r_mic(:,3)-r_src(:,3)',M*N,1,1)]; 
    
    term_sqrt = sqrt( (sum(Mach.*R_vect,2)).^2 + beta2*sum(R_vect.^2,2) ); %intermediate calculation
    
    dr = (-sum(Mach.*R_vect,2) + term_sqrt )/beta2;%distance correction : dr/c = emmision time delay
  
        
    G = exp(Gsign*j*k*dr)./(4*pi*term_sqrt); 
	G=reshape(G,M,N);

%comparaison avec green free field
c=340; %speed of sound in air, m/sbeam
    k=2*pi*f/c;
    
    r_mic=double(r_mic);
    r_src=double(r_src);
    
    R=pdist2(r_mic , r_src);
    
    
    G2 = exp(Gsign*j*k*R)./(4*pi*R);          

%display
figure
imagesc(abs(G))
figure
imagesc(abs(G2))

figure
imagesc(imag(G))
figure
imagesc(imag(G2))
