clear all;
close all;
x=load('erdc11wt.txt') ;
theta=0;
theta1=0;
lat=x(:,1);
lon=x(:,2);
time=x(:,3);
anomaly=x(:,4);
residue=x(:,5);
%total=NBP0207(:,8);
%figure;
count=length(lat);
for i=1:1:count
    if residue(i)>9000
    anomaly(i)=0;
    residue(i)=0;
    end
end
%plot(lat,ano,'r',lat,residue,'b');%,lat,total,'g');
[B,A]= cheby1(3,0.9, 0.1);
 ano=anomaly;% filter(B,A,anomaly);
figure(1);

load synmag1.dat;
aaa= synmag1;
a1=aaa(:,5);

plot(a1);
figure(2);
Y=fft(ano).*exp(1i*pi*theta/180);
Y1=fft(ano).*exp(1i*pi*theta1/180);
%plot(x,real(Y),'b',x,imag(Y),'r');
N=length(ano)
for k=1:N/2-1
    Y(N-k+1)=Y(k+1)';
     Y1(N-k+1)=Y1(k+1)';
end
    
%subplot(4,1,3);
%plot(x,real(Y),'b',x,imag(Y),'r');
xx=ifft(Y);
xx1=ifft(Y1);
%subplot(4,1,4);
%plot(x,real(xx));
plot(lat,real(xx1-mean(xx1)),'r',lat,real(xx-mean(xx)),'b');
sxx=[1.54 2.37];
sa=[578 943];

sas=[1:length(xx)];
sxx1=round(interp1(lat,sas,sxx));
d1=real(xx(sxx1(1):sxx1(2))-mean(xx(sxx1(1):sxx1(2))));

d20=a1(sa(1):sa(2))-mean(a1(sa(1):sa(2)));
xint0=1:length(d20);
xint=linspace(1,length(d20),length(d1));
d2=15*interp1(xint0,d20,xint);

% for i=188:235
%     %d1(186:245)=interp1([186 245],[d1(186) d1(245)],[186:245]);
%     %d1(i)=(d1(i)+40)/8-40;
%     d1(i)=d2(i);
% end
figure(3);
plot(lat(sxx1(1):sxx1(2)),d1,lat(sxx1(1):sxx1(2)),d2,'r');
[Dist,D,k,w,rw,tw]=cdtw(d2,d1,1);

figure(10);

save dd1 d1 d2 ;
for i=1:size(w(:,1))
plot(w(i,1),w(i,2),'.');

hold on;
end
%r=[1 43 88 118 159 189 270];
r=[1  112 152 205 235 342];
rr=size(1,length(r));
for i= 1:length(r)-1
    rr(i)=(w(r(i+1),2)-w(r(i),2))/(w(r(i+1),1)-w(r(i),1));
end