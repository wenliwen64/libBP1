closea 
clear;
figure;
sr = 1; %sampling rate /sec
%t = linspace(0,1,sr);
%x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
load a101;
t = a101(1:3600,1);
x = a101(1:3600,2);
%subplot(3,2,1);                   
[B,A]= cheby1(2,0.9,0.2); %construct the filter
[H,F]=freqz(B,A,1024,sr);
plot(F,abs(H));
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('Chebyshev filter 4th order Frequency response');

%subplot(3,2,2);                  
%[H1,F1]=impz(B,A,1024,sr);
%plot(F1,abs(H1));
%xlabel('Time(s)');
%ylabel('Amplitude');
%title('Chebyshev filter 4th order impulse filter');

%subplot(3,2,3);                  
y = filter(B,A,x);
%plot(t,x,t,y,'r');
%xlabel('Time(s)');
%ylabel('Amplitude');
%title('signals time series');
figure;
%subplot(3,2,4);  
a=fft(x);
b=fft(y);
tt=linspace(0,1,3600);
plot(tt,abs(a),tt,abs(b),'r','markersize',1);
xlim([0 0.5]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('signals in frequency domain');
figure;
plot(t,x,t,y,'r')
%subplot(3,2,5);                  
%zeros=roots(B);
%poles=roots(A);
%zplane(zeros,poles);
%xlabel('Real part');
%ylabel('imaginary');
%title('zeros and poles');



