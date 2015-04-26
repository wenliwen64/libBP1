clear all ;

N1=39;
N2=25;

T0=60;
delay= @(x,y,theta0,c) (x/cosd(theta0)+(y-tand(theta0)*x)*sind(theta0))/c;
c=5e3;%m/s
d1=500;%m
d2=d1;
d=sqrt(d1^2*sind(T0)^2+d2^2*cosd(T0)^2);

load aiii;
load y1;
load y2;
% S=y2;
% I=y1;
t=linspace(0,25,5000);
 S=1*(sin(2*pi*2.5*t)+sin(2*pi*3.5*t-2)+sin(2*pi*4.5*t+2));
 I=100*sin(2*pi*4*t-1);
t0=4;
sr=200;
length=3500;
X=zeros(1,length);
 for i=1:N1
     x(i)=(i*d1-d1)-(N1*d1-d1)/2;
 end
 for j=1:N2
     y(j)=(j*d2-d2)-(N2*d2-d2)/2;
 end
 
    for j=1:N1
        for k=1:N2
        td0=delay(x(j),y(k),T0,c);%AF(i)=AF(i)+AIII(j,k)*exp(1i*2*pi/(c/fc)*(j*d1*cos(theta3(i))+k*d2*sin(theta3(i))));
        td1=delay(x(j),y(k),30,c);
        X=X+AIII(j,k)*(S(round((t0+td0)*sr):round((t0+td0)*sr)+length-1)+I(round((t0+td1)*sr):round((t0+td1)*sr)+length-1));%exp(1i*2*pi*fc*((k-1)*cos(theta3(i))+(j-1)*sin(theta3(i))));
        
        end
    end
    tt=linspace(0,length/sr,length);
    t=linspace(0,25,5000);
    figure(2);
    subplot(3,1,1);
    plot(tt+t0,X,'b');
    xlim([0 25])
    subplot(3,1,2);
    plot(t,S,'b');
    xlim([0 25])
    subplot(3,1,3);
    plot(t,I,'b');
    xlim([0 25]);