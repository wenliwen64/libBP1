clear;
theta=180
x=0:1000;
subplot(4,1,1);
y=sin(x/30);
%z=x*0+exp(i*pi*theta/180);
plot(x,y);
subplot(4,1,2);

Y=fft(y).*exp(i*pi*theta/180);
plot(x,real(Y),'b',x,imag(Y),'r');
N=length(x)
for k=1:N/2-1
    Y(N-k+1)=Y(k+1)';
end
    
subplot(4,1,3);
plot(x,real(Y),'b',x,imag(Y),'r');
xx=ifft(Y);
subplot(4,1,4);
plot(x,real(xx));
