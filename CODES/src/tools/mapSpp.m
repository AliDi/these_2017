function mapSpp(x,y,z)

x=double(x);
y=double(y);
a=double(z);

xi=linspace(min(x),max(x),1000);
yi=linspace(min(y),max(y),1000);
[XI YI]=meshgrid(xi,yi);
ZI = griddata(x,y,z,XI,YI);

figure
contourf(XI,YI,ZI)
xlabel('x')
ylabel('y')

end