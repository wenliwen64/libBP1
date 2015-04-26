clear;
close all;
nel=100;
timesteps=100;
U=zeros(nel,1);
U1=zeros(nel,1);
v=zeros(nel,1);
v1=zeros(nel,1);
a=zeros(nel,1);
a1=zeros(nel,1);
F=zeros(nel,1);
K=ones(nel,1);
d=zeros(nel,timesteps);
for ii=2:nel-1
   if(ii>50)
    K(ii)=0.3;
   end
    
end
dt=1.0*(1/1)^0.5;
Damp=0.00;

figure;
xx=linspace(0,1,nel);
M=1;
U(25)=1;
for kk=1:timesteps
for ii=2:nel-1
    F(ii)=K(ii)*(U(ii+1)-U(ii))-K(ii)*(U(ii)-U(ii-1));
end
F(1)=K(1)*(U(2)-U(1));
F(nel)=-K(nel)*(U(nel)-U(nel-1));
a1=(F-Damp*v)/M;
v1=(a+a1)/2*dt+v;
U1=U+v1*dt+dt^2/2*a1;
a=a1;
v=v1;
U=U1;
d(:,kk)=U;
% if mod(kk-1,1000)==0
%         subplot(6,4,(kk-1)/1000+1);
%         plot(xx,U,'r');
%         ylim([-0.5 0.5]);
% end
end

lat=zeros(nel,timesteps);
lon=zeros(nel,timesteps);
for iii=1:nel
        lat(iii,:)=(1:timesteps)*dt;
end

for iii=1:timesteps
        lon(:,iii)=xx;
end
      figure(1);
        H=pcolor(lon, lat, d);

%colormap(gray);

        shading interp;
        colorbar;       