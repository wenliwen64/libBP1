close all;
clear all;
scale=500;
ranscal=25;
b0=linspace(0,1,ranscal);
b=linspace(0,1,scale);
x0=rand(ranscal);

for i=1:ranscal
    lon1(i,:)=b0;
end
for i=1:ranscal
    lat1(:,i)=b0;
end


for i=1:scale
    lon2(i,:)=b;
end
for i=1:scale
    lat2(:,i)=b;
end
x=interp2(lon1,lat1,x0,lon2,lat2);
xx=zeros(2*scale);
for i=1:scale
    for j=1:scale
    xx(scale+i,scale+j)=x(i,j);
    xx(scale+1-i,scale+j)=x(i,j);
    xx(scale+i,scale+1-j)=x(i,j);
    xx(scale+1-i,scale+1-j)=x(i,j);
    end
end
for i=1:2*scale
    lon(i,:)=linspace(-1,1,scale*2)*100;
end
for i=1:2*scale
    lat(:,i)=linspace(-1,1,scale*2)*100;
end
figure(1)
pcolor(lon,lat,xx);
shading interp;
colorbar;
yy=fft2(xx);
for i=1:2*scale
    lon(i,:)=linspace(0,10*pi,scale*2);
end
for i=1:2*scale
    lat(:,i)=linspace(0,10*pi,scale*2);
end
figure(2)
pcolor(lon,lat,imag(yy));
colorbar;
shading interp;
figure(3)
pcolor(lon,lat,real(yy));
colorbar;
shading interp;
caxis([0 1000]);