clear all;
close all;
m=5000;
time=2;
ttt=linspace(0,time,m*time)*2*pi;
kkk=linspace(0,m,m*time);
%tt=rand(m,1);
t=cos(200*ttt);%+cos(2*ttt);
% for i=1:m
%     t(i+m)=tt(i);
%     t(m+1-i)=tt(i);
% end
%tt=linspace(-2*pi,2*pi,719);
%CFUN= @(b) cos(b);
%xx=fft(CFUN(t));
%yy=ifft(xx);

figure(1);
subplot(2,1,1);
plot(ttt,t);
subplot(2,1,2);
plot(kkk,real(fft(t)),'r',kkk,imag(fft(t)),'b');
xlim([0 1000]);
%figure(2);
%plot(tt,xcorr(CFUN(t)),'b',t,ifft(xx.*xx),'r');
