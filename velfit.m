close all
clear all
P=[7.29320740145629e-60,-2.25713355056818e-54,3.15780784770224e-49,-2.64134885422995e-44,1.47242246990781e-39,-5.77410811052121e-35,1.63959228354946e-30,-3.42090336118034e-26,5.27212978329383e-22,-5.99116415754033e-18,4.98013313198646e-14,-2.98197497509444e-10,1.25141562627232e-06,-0.00349386067308002,5.89649663355069,1273.93757459125;]

%parkersfield velocity model

CFUN= @(b) P(16)+P(15)*b+P(14)*b.^2+P(13)*b.^3+P(12)*b.^4+P(11)*b.^5+(P(10)+P(9)*b+P(8)*b.^2+P(7)*b.^3+P(6)*b.^4+P(5)*b.^5+P(4)*b.^6+P(3)*b.^7+P(2)*b.^8+P(1)*b.^9).*b.^6;
% clear all;
% 
% close all;
% load b0;
a=0:0.1:30;
b=CFUN(a*1000);
b(200:301)=b(200:301)+(0:101)*3.5;
figure(1);
plot(b,-a,'.');
P=polyfit(a,b,15);
a1=0:0.1:30;
b1=polyval(P,a1);
figure(2);
plot(b,-a1);
