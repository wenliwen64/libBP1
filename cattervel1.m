close all;
clear all;
scale=401;
length=10;
nyq=pi/(length/scale);
ranscal=401;
a=1;
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
% xx=zeros(2*scale);
% for i=1:scale
%     for j=1:scale
%     xx(scale+i,scale+j)=x(i,j);
%     xx(scale+1-i,scale+j)=x(i,j);
%     xx(scale+i,scale+1-j)=x(i,j);
%     xx(scale+1-i,scale+1-j)=x(i,j);
%     end
% end
for i=1:scale
    lon0(i,:)=linspace(0,1,scale)*length;
end
for i=1:scale
    lat0(:,i)=linspace(0,1,scale)*length;
end
% figure(1)
% pcolor(lon0,lat0,x);
% shading interp;
% colorbar;
temp=zeros(scale);
yy=zeros(scale);
% for i=1:scale9
%     temp(i,:)=fft(x(i,:));
% end
% for i=1:scale
%     yy(:,i)=fft(temp(:,i));
% end
yy=fft2(x);













yyy=zeros(scale);
yyy(1,:)=yy(1,:);
yyy(:,1)=yy(:,1);


for i=1:scale
    lon(i,:)=linspace(0,2*nyq,scale);
end
for i=1:scale
    lat(:,i)=linspace(0,2*nyq,scale);
end
for i=2:scale
    for j=2:(scale+1)/2
        if i<=(scale+1)/2
            c=2;
            b=2;
        end
        if i>(scale+1)/2
            c=scale;
            b=2;
        end
        
            
        %theta=atan((imag(yy(i,j)))/(real(yy(i,j))));
        kr=((i-c)^2+(j-b)^2)^0.5*2*nyq/scale;
        %yyy(i,j)=a^2/(1+kr^2*a^2)*(cos(theta)+sin(theta)*(-1)^0.5);
        %yyy(i,j)=((imag(yy(i,j)))^2+(real(yy(i,j)))^2)^0.5*(cos(theta)+sin(theta)*1i);    
        yyy(i,j)=yy(i,j)/abs(yy(i,j))*(1e5*a^2/(1+kr^2*a^2*1e-2))^0.5;
                %yyy(i,j)=yy(i,j)/abs(yy(i,j))*(a^2/(1+kr^2*a^2))^0.5;
        %yyy(i,j)=yy(i,j);
    end
end
% for i=1:scale
%      if i<=(scale+1)/2
%             c=1;
%             b=1;
%         end
%         if i>(scale+1)/2
%             c=scale;
%             b=1;
%         end
%    
%      kr=((i-c)^2+(1-b)^2)^0.5*2*nyq/scale;
%         %yyy(i,j)=a^2/(1+kr^2*a^2)*(cos(theta)+sin(theta)*(-1)^0.5);
%         %yyy(i,j)=((imag(yy(i,j)))^2+(real(yy(i,j)))^2)^0.5*(cos(theta)+sin(theta)*1i);    
%         yyy(i,1)=yy(i,1)/abs(yy(i,1))*(1e5*a^2/(1+kr^2*a^2*1e-2))^0.5;
%         
%      
% end
% for i=1:scale
%      if i<=(scale+1)/2
%             c=1;
%             b=1;
%         end
%         if i>(scale+1)/2
%             c=1;
%             b=scale;
%         end
%    
%      kr=((1-c)^2+(i-b)^2)^0.5*2*nyq/scale;
%         %yyy(i,j)=a^2/(1+kr^2*a^2)*(cos(theta)+sin(theta)*(-1)^0.5);
%         %yyy(i,j)=((imag(yy(i,j)))^2+(real(yy(i,j)))^2)^0.5*(cos(theta)+sin(theta)*1i);    
%         yyy(1,i)=yy(1,i)/abs(yy(1,i))*(1e5*a^2/(1+kr^2*a^2*1e-2))^0.5;
%         
%      
% end
for i=2:scale
    for j=2:(scale+1)/2
        yyy(scale-i+2,scale-j+2)=yyy(i,j)';
    end
end


% 
% for i=2:(scale+1)/2untitled
%     for j=2:(scale+1)/2
%      yyy(scale-i+2,j)=yy(i,j)';
%      yyy(i,scale+2-j)=yy(i,j)';
%      yyy(scale+2-i,scale+2-j)=yy(i,j);
%      yyy(i,j)=yy(i,j);
%     end
% end
xx=ifft2(yyy);
figure(2)
subplot(2,3,1);
pcolor(lon,lat,real(yy));
colorbar;
shading interp;
caxis([0 1000]);


subplot(2,3,2)
pcolor(lon,lat,imag(yy));
colorbar;
shading interp;
%caxis([0 1000]);

subplot(2,3,3)
pcolor(lon,lat,abs(yy));
colorbar;
shading interp;
caxis([0 1000]);

subplot(2,3,4)
pcolor(lon,lat,real(yyy));
colorbar;
caxis([0 0.01]);
shading interp;



subplot(2,3,5)
pcolor(lon,lat,imag(yyy));
colorbar;
shading interp;
%caxis([-1 1]);4

subplot(2,3,6)
pcolor(lon,lat,abs(yyy));
colorbar;
shading interp;
caxis([0 0.01]);

figure(3)
subplot(2,1,1)
pcolor(lon0,lat0,real(x));
colorbar;
shading interp;
title('random media');

% subplot(2,2,2)
% pcolor(lon0,lat0,imag(x));
% colorbar;
% shading interp;


subplot(2,1,2)
pcolor(lon0,lat0,real(xx));
colorbar;
shading interp;
caxis([0 1 ])
title('Von Karman media(2d FFT and filtered from the first figure)');
% subplot(2,2,4)
% pcolor(lon0,lat0,imag(xx));
% colorbar;
% shading interp;
% caxis([0 1 ])



figure(6);
y1=yy(:,2);
for i=2:scale
    
    test(i)=-imag(y1(scale+2-i));
end
    test(1)=-imag(y1(1));
plot([1 :scale],imag(y1),[1 :scale], test);

figure(7);
y1=yy(:,2);
for i=2:scale
    
    test(i)=real(y1(scale+2-i));
end
    test(1)=real(y1(1));
plot([1 :scale],real(y1),[1 :scale], test);




figure(8);
krr=linspace(0,2*nyq/scale*200*1.414,100);
plot(krr,1e5*a^2./(1+krr.^2*a^2*1e-2));

