
function ss=YoffeAnalytic(tt,tr,ts,Dmax)
%tt is the input time
%tr is the rise time
%ts is the process zone time
%Dmax is just a coefficient of the max amplitude
%ss is the output yoffe function
N=length(tt);
ss=zeros(1,N);
for i=1:N
    t=tt(i);
C1=(0.5*t+0.25*tr)*sqrt(t*(tr-t))+(t*tr-tr*tr)*asin(sqrt(t/tr))-0.75*tr^2*atan(sqrt((tr-t)/t));
C2=3/8*pi*tr^2;
C3=(ts-t-0.5*tr)*sqrt((t-ts)*(tr-t+ts))+tr*(2*tr-2*t+2*ts)*asin(sqrt((t-ts)/tr))+3/2*tr^2*atan(sqrt((tr-t+ts)/(t-ts)));
C4=(-ts+0.5*t+0.25*tr)*sqrt((t-2*ts)*(tr-t+2*ts))+tr*(-tr+t-2*ts)*asin(sqrt((t-2*ts)/tr))-0.75*tr^2*atan(sqrt((tr-t+2*ts)/(t-2*ts)));
C5=pi/2*tr*(t-tr);
C6=pi/2*tr*(2*ts-t+tr);
K=2/(pi*tr*ts*ts);
if t<0% assume tr>2*ts
    s=0;
elseif t<ts
    s=C1+C2;
elseif t<2*ts
    s=C1-C2+C3;
elseif t<tr
    s=C1+C3+C4;
elseif t<tr+ts
    s=C5+C3+C4;
elseif t<tr+2*ts
    s=C4+C6;
 else
    s=0;
end
ss(i)=Dmax*s*K;
end
ss=ss/sum(ss);
end
% function test()
% figure(1);
% t=-5:0.1:15;
% tr=5;
% ts=1;
% Dmax=1;
% plot(t,Yoffe1(t,tr,ts,Dmax));
% end