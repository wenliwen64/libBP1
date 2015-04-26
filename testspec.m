close all;
clear all;
load xx2;
x=diff(diff(xx2(2,:)));
%x=xx2(2,:);
sr=5;
t=(0:length(x)-1)/sr;
figure(1);
plot(t,x);
%[B,A]=butter(4,[0.5 1.5]/(sr/2));
%x1=filter(B,A,x);
x1=x;
%xlim off;
xlim([100 300]);
t0=140;
t1=160;
t2=180;
N=(t1-t0)*sr+1;
ff=linspace(0,sr-1/(t1-t0),N);
fs=abs(fft(x(t0*sr:t1*sr)));
fs1=abs(fft(x1(t0*sr:t1*sr)));
figure(3)
%cs=abs(cmtm(x(t0*sr:t1*sr),x(t0*sr:t1*sr),2,3,1:(t1-t0)*sr+1,1));
cs=abs(cmtm2(x(t0*sr:t1*sr),x(t0*sr:t1*sr),2));
[ps w]=pmtm(x(t0*sr:t1*sr),2);
ps=abs(ps);
w=w*sr/2/pi;
figure(2);
plot(ff,fs/max(fs),'r',ff,fs1/max(fs1),'b',ff,cs/max(cs),'g',w,ps/max(ps),'r');
xlim([0 sr/2]);


