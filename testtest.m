t=linspace(0,2*pi,100);
x=sin(t);
y=ifft(fft(x)*exp(1i*pi/4));
plot(t,x,'b',t,y,'r');