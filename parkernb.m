clear all;
close all;
% least mean square adaptive array 
k=@(lamda,theta) 2*pi/lamda*[cos(theta) sin(theta)];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
delay= @(x,y,theta0,c) (x/cos(theta0)+(y-tan(theta0)*x)*sin(theta0))/c;
d=0.25;%km
nel=14;
vp=5;%km
r=zeros(nel,2);%coordinates of the array element
%r(:,1)=[-371 -298 -192 38 -23 0 23 177 213 260 98 -72 -131 -195]/1000;
%r(:,2)=[-302 -208 -284 -180 -13 0 -9 98 235 411 213 328 393 356]/1000;
r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
time=25;%total time in seconds
samplerate=200;
length=4000;
t=linspace(0,time,samplerate*time);
L0=5;
L=L0+5;
f0=10;%Hz
f1=6;
f2=6;
SNR=1000000;
theta0=110*pi/180;
theta1=130*pi/180;
theta2=60*pi/180;
A=1;
% S=zeros(1,samplerate*time);
% for i=1:samplerate *time
%     if t>=0
% S(i)=A*(sin(2*pi*f0*t(i))+sin(2*pi*1.2*t(i)-2)+sin(2*pi*0.8*t(i)+2));
%     end
% end
S=A*(sin(2*pi*f0*t)+sin(2*pi*8*t-2)+sin(2*pi*7*t+2));
I1=30*sin(2*pi*f1*t-1);
I2=30*(sin(2*pi*f2*t-3));
load y1;
load y2;
% S=y2;
% I1=y1;
% I2=y1;
t0=2;

N=zeros(nel,length);
X=zeros(nel,length);
for i=1:nel
    N(i,:)=0;%1/SNR*randn(1,length);
    td0=delay(r(i,1),r(i,2),theta0,vp);%AF(i)=AF(i)+AIII(j,k)*exp(1i*2*pi/(c/fc)*(j*d1*cos(theta3(i))+k*d2*sin(theta3(i))));
     td1=delay(r(i,1),r(i,2),theta1,vp);
     td2=delay(r(i,1),r(i,2),theta2,vp);
     %X=v(r,k(vp/f0,theta0))*S+N+v(r,k(vp/f1,theta1))*I1+v(r,k(vp/f2,theta2))*I2;
     X(i,:)=S(round((t0+td0)*samplerate):round((t0+td0)*samplerate)+length-1)+I1(round((t0+td1)*samplerate):round((t0+td1)*samplerate)+length-1)+I1(round((t0+td2)*samplerate):round((t0+td2)*samplerate)+length-1);
end
%X=v(r,k(vp/f0,theta0))*S+N+v(r,k(vp/f1,theta1))*I1+v(r,k(vp/f2,theta2))*I2;


Rxx=zeros(nel,nel);
for i=1:length
    Rxx=Rxx+X(:,i)*X(:,i)';
end
Wopt=A^2*inv(Rxx)*(v(r,k(vp/21.3,theta0)));


%plot(t,X(1,:))
Y=Wopt'*X;
Y1=ones(1,nel)*X;
%e=Y-S;

%MSE=mean(abs(e));
 tt=linspace(t0,length/samplerate+t0,length);
res=1000;
theta3=linspace(0,-pi,res);
AF=zeros(1,res);
for i=1:res
    vv=v(r,k(vp/f0,theta3(i)));
    for j=1:nel
    AF(i)=AF(i)+Wopt(j)*vv(j);
    end
end
figure(3);
N=2;
% subplot(N,1,1)
% plot(theta3*180/pi+180,10*log10(abs(AF).^2/nel^2));

subplot(N,1,1)
plot(t,S/max(S),'b',tt,Y/max(Y),'r');
xlim([L0 L]);
% subplot(N,1,3)
% plot(tt,Y/(max(Y)),'r',tt,Y1/max(Y1),'g');
% xlim([L0 L]);
subplot(N,1,2)
plot(t,I1/(max(I1)),'black',tt,Y1/max(Y1),'g');
xlim([L0 L]);


step=0;
Wadp=ones(1,nel);
Yadp=zeros(1,samplerate*time);
%%%%%%%% LMS adptive
% for i=1:length
%     
%     Rxx=X(:,i)*X(:,i)';
%     step=0.002/max(eig(Rxx));
%     Wadp=Wadp+step*(A^2*v(r,k(vp/f0,theta0))-Rxx*Wadp')';
%     Yadp(i)=Wadp*X(:,i);
% end
% subplot(N,1,4)
% plot(tt,Yadp/(max(Yadp)),'r',tt,Y1/max(Y1),'g');
% xlim([0 30]);
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

