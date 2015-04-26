clear;
load loc_dat_1;
lat=loc_dat_1(:,1);
lon=loc_dat_1(:,2);
z=loc_dat_1(:,3)/1000;
x=loc_dat_1(:,4);
y=loc_dat_1(:,5);
t=loc_dat_1(:,6);

x0=0;
y0=0;
z0=-10;
t0=0;
v=6;
N=length(lat);
for kk=1:1:10
for k=1:1:N
ti0(k)=t0+sqrt((x(k)-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v;
%pdx(k)=(sqrt((x(k)+0.1-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v-sqrt((x(k)-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v)/0.1;
pdx(k)=-(x(k)-x0)/sqrt((x(k)-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v;
pdy(k)=-(y(k)-y0)/sqrt((x(k)-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v;
pdz(k)=-(z(k)-z0)/sqrt((x(k)-x0)^2+(y(k)-y0)^2+(z(k)-z0)^2)/v;
one(k)=1;
dt(k)=t(k)-ti0(k);
end
A=[pdx;pdy;pdz;one]';
C=inv(A'*A);
m=C*A'*dt';
x0=x0+m(1);
y0=y0+m(2);
z0=z0+m(3);
t0=t0+m(4);
delta_m1=0;
delta_m2=0;
delta_m3=0;
delta_m4=0;
for kkk=1:1:N
delta_m1=delta_m1+C(1,1)*(t(kkk)-ti0(kkk))^2/(N-4);
delta_m2=delta_m2+C(2,2)*(t(kkk)-ti0(kkk))^2/(N-4);
delta_m3=delta_m3+C(3,3)*(t(kkk)-ti0(kkk))^2/(N-4);
delta_m4=delta_m4+C(4,4)*(t(kkk)-ti0(kkk))^2/(N-4);
end
delta_m1=delta_m1^0.5;
delta_m2=delta_m2^0.5;
delta_m3=delta_m3^0.5;
delta_m4=delta_m4^0.5;
[x0 y0 z0 t0]
[delta_m1 delta_m2 delta_m3 delta_m4]
end
