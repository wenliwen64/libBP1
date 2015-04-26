t=linspace(-2*pi,2*pi,100);
tt=linspace(-2*pi,0,50);
S=sin(t);
vp=1;
f=1;
lamda=vp/f;
t0=0;
S1=S*cos(2*pi*t0*1/f0);
S2=S(1+t0*50:1+t0*50+50-1);
figure(1);
plot(t,S1,'b',tt,S2,'r',t,S,'g')



