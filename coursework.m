clear;
figure;

sr1 = 1 %sampling rate /sec
sr20 = 20
sr100 = 100
d = load('L25.m', 'v1');
xt=d(1:86400*sr1,1);
yt=d(1:86400*sr1,2);
t = linspace(0,sr1,3600*sr1);
%x=
z1=1.1*exp(i*2*pi*0.16/sr1);
cz1=conj(z1);
z0=1.35*exp(i*2*pi*0.16/sr1);
cz0=conj(z0);
B=[z0*cz0 -(z0+cz0) 1];
A=[z1*cz1 -(z1+cz1) 1];
[H1,F1]=freqz(B,A,1024,sr1);



subplot(2,2,1)
plot(F1,abs(H1));
xlim([0 0.5])
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('filter frequency response');

%y=filter(B,A,x);
subplot(2,2,2)
zplane(roots(B),roots(A));
title('postions of poles and zeros on the Z-plane');

%figure;

subplot(2,2,[3 4]);  
fftyt=abs(fft(yt(1:3600*sr1)));
yf1=filter(B,A,yt);
yf=yf1-yt;
fftyf=abs(fft(yf(1:3600*sr1)));

plot(t(2:3600*sr1),fftyt(2:3600*sr1),'b',t(2:3600*sr1),fftyf(2:3600*sr1),'r');
xlim([0 0.5]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('signals in frequency domain');

% subplot(5,2,[5 6])
% plot(xt,yt);
% xlim([10 1000]);
% xlabel('Time(s)');
% ylabel('Amplitude');
% title('signals in time domain');
% 
% 
% 
% 
% subplot(5,2,[ 7 8 ])
% plot(xt,yf)
% xlim([0 1000]);
% xlabel('Time(s)');
% ylabel('Amplitude');
% title('signals at 1/6 Hz frequency in time domain ');
% 
% 
% f(1:24)=0
% for k=1:1:24
%     for kk=1:1:3600*sr1
%         f(k)=f(k)+yf((k-1)*3600*sr1+kk)^2;
%     end
% end
% %figure;
% subplot(5,2,[ 9 10 ])
% plot([0:23],f,'r*-')
% xlabel('Time(hour)');
% ylabel('Amplitude');
% title('sum of the energy of the signals in an hour');
% 
