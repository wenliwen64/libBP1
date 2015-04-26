clear all;
close all;
% least mean square adaptive array 
k=@(lamda,theta) 2*pi/lamda*[cos(theta) sin(theta)];% wave vector
v=@(r,k) exp(1i*(r*k'));%steering vector
delay= @(x,y,the,c) (x/cos(the)+(y-tan(the)*x)*sin(the))/c;
d=0.125;%km
nel=9;
vp=5;%km
r=zeros(nel,2);%coordinates of the array element
r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
time=25;%total time in seconds
samplerate=500;
length=4000;
length2=3000;
t0=2;
t=linspace(0,time,time*samplerate);
 %t=linspace(t0,length/samplerate+t0,length);
L0=0;
L=1;
fm0=1;
f0=20;%Hz
f1=20.0;
f2=20.0;
SNR=10000;
theta0=90*pi/180;
theta1=160*pi/180;
theta2=140*pi/180;
A=1;
% S=zeros(1,samplerate*time);
% for i=1:samplerate *time
%     if t>=0
% S(i)=A*(sin(2*pi*f0*t(i))+sin(2*pi*1.2*t(i)-2)+sin(2*pi*0.8*t(i)+2));
%     end
% end
S=hilbert(sin(2*pi*t*fm0).*sin(2*pi*f0*t));
I1=100*hilbert(cos(2*pi*f1*t-3));
I2=100*hilbert((sin(2*pi*f2*t-3)));
load y1;
load y2;
% S=y2;
% I1=10*y1;
% I2=10*y1;


N=zeros(nel,time*samplerate);
X=zeros(nel,time*samplerate);
% XX=zeros(nel,samplerate*time);
% td0=zeros(1,nel);
% td1=zeros(1,nel);
% td2=zeros(1,nel);
for i=1:nel
     N(i,:)=1/SNR*randn(1,time*samplerate);
%     td0(i)=delay(r(i,1),r(i,2),theta0,vp);%AF(i)=AF(i)+AIII(j,k)*exp(1i*2*pi/(c/fc)*(j*d1*cos(theta3(i))+k*d2*sin(theta3(i))));
%      td1(i)=delay(r(i,1),r(i,2),theta1,vp);
%      td2(i)=delay(r(i,1),r(i,2),theta2,vp);
%      %X=v(r,k(vp/f0,theta0))*S+N+v(r,k(vp/f1,theta1))*I1+v(r,k(vp/f2,theta2))*I2;
%      XX(i,:)=S(round((t0+td0(i))*samplerate):round((t0+td0(i))*samplerate)+length-1)+I1(round((t0+td1(i))*samplerate):round((t0+td1(i))*samplerate)+length-1)+I2(round((t0+td2(i))*samplerate):round((t0+td2(i))*samplerate)+length-1);
 end
X=v(r,k(vp/f0,theta0))*S+N+v(r,k(vp/f1,theta1))*I1+v(r,k(vp/f2,theta2))*I2;


Rxx=zeros(nel,nel);
for i=1:time*samplerate
    Rxx=Rxx+X(:,i)*X(:,i)';
end

%     YY=inv(Rxx)*X;
%     YYY=zeros(1,length2);
% for i=1:nel
%     l,samplerate*time);
% td0=zeros(1,nel);
% td1=zeros(1,nel
%     
%     YYY=YYY+YY(i,round((t0-td0(i))*samplerate):round((t0-td0(i))*samplerate)+length2-1);
%     
%     
%     
% end






Wopt=A^2*inv(Rxx)*(v(r,k(vp/f0,theta0)));


%plot(t,X(1,:))
Y=(Wopt'*X);
Y1=(ones(1,nel)*X);
%e=Y-S;

%MSE=mean(abs(e));
%  tt=linspace(t0,length/samplerate+t0,length);
%  ttt=linspace(2*t0,length2/samplerate+2*t0,length2);
res=1000;
theta3=linspace(0,pi,res);
AF=zeros(1,res);
for i=1:res
    vv=v(r,k(vp/f0,theta3(i)));
    for j=1:nel
    AF(i)=AF(i)+Wopt(j)*vv(j);
    end
end
figure(3);
N=4;
subplot(N,1,1);
plot((pi-theta3)*180/pi,10*log10(abs(AF).^2/nel^2));

subplot(N,1,2);
plot(t,S/max(S),'b',t,Y/(max(Y)),'r',t,Y1/(max(Y1)),'g');
xlim([L0 L]);
subplot(N,1,3)
plot(t,S/max(S),'b',t,I1/max(I1),'r');
xlim([L0 L]);
subplot(N,1,4)
plot(t,I1/max(I1),'b',t,Y/(max(Y)),'r');
% step=0;
% Wadp=ones(1,nel);
% Yadp=zeros(1,length);
% %%%%%%%% LMS adptive
%plot(tt,X(9,:),'r',t,XX(9,:),'b');ttt,YYY/(max(YYY)),'r'
xlim([L0 L]);
% step=0;
% Wadp=ones(1,nel);
% Yadp=zeros(1,length);
% %%%%%%%% LMS adptive
% for i=1:length
%     
%     Rxx=X(:,i)*X(:,i)';
%     step=0.002/max(eig(Rxx));
%     Wadp=Wadp+step*(A^2*v(r,k(vp/f0,theta0))-Rxx*Wadp')';
%     Yadp(i)=Wadp*X(:,i);
% end
% subplot(N,1,4)
% plot(t,Yadp/(max(Yadp)),'r',tt,Y1/max(Y1),'g');
% xlim([L0 L]);
% subplot(N,1,5)
% plot(tt,Yadp/(max(Yadp)),'r',tt,Y1/max(Y1),'g');
% xlim([L0 L]);
% subplot(N,1,4)
% plot(t,-Yadp/(max(Yadp)),'r',t,Y1/max(Y1),'g');
% xlim([100 130]);
% subplot(N,1,5)
% plot(t,-Yadp/(max(Yadp)),'r',t,Y1/max(Y1),'g');
% xlim([370 400]);
% subplot(N,1,8)
% plot(t,-Yadp/(max(Yadp)),'r',t,Y1/max(Y1),'g');
% xlim([120 150]);

