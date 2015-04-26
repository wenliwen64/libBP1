G=6.67E-11
rou=4500
Cp=1e7*1e-7*1000
R=linspace(0,6370e3,100)
T=4*G*rou*pi*R.^2/5/Cp+300;
plot(R/1000,T);
xlabel('radius (km)');
ylabel('temperature (K)');
title('profile of the proto earth temperature (conductivity = inf)');
