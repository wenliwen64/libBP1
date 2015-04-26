figure;
sr = 400 %sampling rate /sec
t = linspace(0,1,sr);
x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
z1=1.02*exp(i*2*pi*55/sr);
cz1=conj(z1);
z0=1.0005*exp(i*2*pi*55/sr);
cz0=conj(z0);
%zplane(z0,z1);
B=[z0*cz0 -(z0+cz0) 1];
A=[z1*cz1 -(z1+cz1) 1];
[H1,F1]=freqz(B,A,1024,sr);
[H2,F2]=impz(B,A,1024,sr);
%plot(F2,abs(H2));
y=filter(B,A,x);
%plot(t,fft(x),t,fft(y1),'r')
%xlim([0 0.5]);



subplot(3,2,1);                   
[H1,F1]=freqz(B,A,1024,sr);
plot(F1,abs(H1));
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('Notch filter Frequency response');


subplot(3,2,2);                  
[H2,F2]=impz(B,A,1024,sr);
plot(F2,abs(H2));
xlabel('Time(s)');
ylabel('Amplitude');
title('Notch filter impulse response');


subplot(3,2,3);                  
y = filter(B,A,x);
plot(t,x,t,y,'r');
xlabel('Time(s)');
ylabel('Amplitude');
title('signals time series');

subplot(3,2,4);                  
plot(t,fft(x),t,fft(y),'r','markersize',2);
xlim([0 0.5]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('signals in frequency domain');

subplot(3,2,5);                  
xlabel('Real part');
ylabel('imaginary');
title('zeros and poles');
zplane(roots(B),roots(A));

