figure;
title('orbit4096')
N=4096
p=orbit(:,1);
a1=orbit(:,2);
b1=orbit(:,3);
dt=p(2)-p(1)
samplerate=1/dt
nyq=samplerate/2
x=linspace(-nyq,nyq,N);
    

%plot(orbit(:,2),orbit(:,3))
%a1+b1*i
subplot(3,1,1);
semilogy(x,abs(fftshift(a1+b1*i)));
xlabel('Frequency(1/earth year)')  ;
title('Straight Fourier Transform');





%xlim([-5,0])
subplot(3,1,2)
x1=linspace(-nyq,nyq,2*N-1);
y1=abs(fft(xcorr(a1+b1*i,a1+b1*i),2*N-1))
x1=x1
semilogy(x1,y1);
xlabel('Frequency(1/earth year)')  ;
title('Autocorrelation');








subplot(3,1,3);
y2=fft(arburg(a1+b1*i,3200))
for j=1:1:3201
    temp(j)=abs(y2(j))^2;%(j)*conj(y2(j));
    pwr1(j)=1/temp(j);
end
%pwr=1/temp;

%power1=1/(y2(1:801).^2);

x2=linspace(-nyq,nyq,3201)
semilogy(x2,pwr1);
%y2=fft(ar(a1+b1*i,400,'burg'))
xlabel('Frequency(1/earth year)')  ;
title('AR model');
