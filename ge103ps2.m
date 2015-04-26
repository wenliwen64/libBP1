F0=1367
sigma=5.67e-8
tou=0.63

T=linspace(200,350,1000);
Tt=F0/4*(1-0.45+0.25*tanh((T-272)/23))-sigma*T.^4/(1+tou);
plot(T,Tt,'b',T,T*0,'r');
title('solutions for T');
xlabel('T(K)');
ylabel('CdT/dt');
for p=1:1:999
Ttt(p)=(Tt(p+1)-Tt(p))/(T(p+1)-T(p));
end
hold on;
plot(T(1:999),Ttt,'b');
T0=(F0*(1+tou)/(4*sigma))^(1/4);
for p=1:1:1000
T00(p)=T(p)*(1-0.45+0.25*tanh((T(p)-272)/23))^(-1/4);
end;
figure;
plot(T,T00,'b',T,T*0+T0,'r');
title('T0 vs T');
xlabel('T');
ylabel('T0');
Fn=313.8^4*4*sigma/(1+0.63);
dFn=(Fn-F0)/F0;
Fb=322.1^4*4*sigma/(1+0.63);
dFb=(Fb-F0)/F0;
T1=(F0*(1-0.3)/(4*sigma))^(1/4);
touT=0.5+0.13*exp(5413*(1/288-T.^(-1)))
T11=T.*(1+touT).^(-1/4);
figure;
plot(T,T11,'b',T,T*0+T1,'r');
title('T1 vs T');
xlabel('T');
ylabel('T1');
Fr=261^4*4*sigma/(1-0.3);
dFr=(Fr-F0)/F0;
Fv=F0/0.732^2;
Te=(F0/sigma/pi)^(1/4);


















