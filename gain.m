nel=6;
vp=5;%km/s
kk=1;
w=[ 4.72 10.92 15.97 15.97 10.92 4.72]/15.97;
for d=[1 2.5 5] ;%km
figure(1);
f0=3;%Hz
res=1000;
theta=linspace(-pi/2,pi/2,res);
tou=d*sin(theta)/vp;
lamda0=vp/f0;
phi=2*pi*d/lamda0*sin(theta);
i=1:res;
A=zeros(1,res);
B=zeros(1,res);
for i=1:res
    for j=1:nel
    A(i)=A(i)+exp(1i*((j-1)*phi(i)));
    end
end
for i=1:res
    for j=1:nel
    B(i)=B(i)+w(j)*exp(1i*((j-1)*phi(i)));
    end
end
subplot(3,1,kk);
kk=kk+1;
plot(theta*180/pi,10*log10(abs(A).^2/nel^2),'b',theta*180/pi,10*log10(abs(B).^2/nel^2),'r');
ylabel('signal gain (dB)');
xlabel(['sigal angle from the broadside ' 'd=' num2str(d) 'km (interelement distance)' ]);
ylim([-60 inf]);
end