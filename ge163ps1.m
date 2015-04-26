rou=3300;
Cp=1000;
k=3;
T1T0=1539;
T0=273;
v=0.075/365/24/3600;
L=100000;
x=(1:200)*1e6*365*24*3600*v;
R=rou*Cp*v*L/2/k
sum=1:200;
for m=1:200
sum(m)=0;
for n=-100:1:100
sum(m)=sum(m)+2*exp(-((n^2*pi^2+R^2)^0.5-R)*x(m)/L);
%sum1(n+101)=(2*exp(-((n^2*pi^2+R^2)^0.5-R)*x/L));
end
qstz(m)=3*((T1T0/L)*(sum(m)-2)+(T1T0/L));
end
plot(x/(1e6*365*24*3600*v),qstz);
title('qstz vs age at L=100000m ');
xlabel('age (Ma)');
ylabel('qs(t)z (W/m2)');


