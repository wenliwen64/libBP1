clear;
Tz=[140.85 65.19 62.16 54.13 53.98 51.28 50.95 50.40 49.94 49.47]
Sy=[0.84 0.99 0.98 1.13 1.11 1.31 1.23 1.36 1.38 1.28]
beta=[1 1 1 1 ]
beta=nlinfit(Tz,Sy,'hyperbola',beta)
h=linspace(49,140,100)
%plot(Tz,Sy,'b.',h,hyperbola(beta,h),'r-')
xlabel('Ti/Zr');
ylabel('Sm/Yb');
title('best-fitting hyperbola')
%figure;
h1=linspace(-1,1,100)
for k=1:1:length(h1)
    v(k)=1/h1(k);
    v1(k)=100/h1(k);
end
%plot(h1,v,'r',h1,v1)

t=1:1:30000;

Mgd(1)=1e-7;
for k=1:1:29999
    Mgd(k+1)=1/1.8e18*3.6e15*1e-8*(2.43+1.5*sin(-0.1*pi*t(k)/5e3))+(1-1/1.8e18*3.6e14)*Mgd(k);
end

Mgd1(1)=1e-7;
for k=1:1:29999
    Mgd1(k+1)=1/1.8e18*3.6e15*1e-8*(2+0*sin(pi/2-0.2*pi*t(k)/5e3))+(1-1/1.8e18*3.6e14)*Mgd1(k);
end

Mgd2(1:2000)=1e-7;
Mgd2(2001:9000)=Mgd1(1:7000);
Mgd2(9001:30000)=Mgd(9001:30000);



%figure;
plot(t,Mgd2*1000)
xlabel('time(yr)');
ylabel('Mg/Ca');
title('Mg/Ca ratio in the sediments as a function of time(inverse calculation) ')