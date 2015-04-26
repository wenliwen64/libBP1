

clear all;
close all;
%case A 
nel=20;
P=8;
omiga=2*3.14*2.5;
h=[ones(1,nel)*0.5];








kzz=zeros(P+1,P+1,nel);
mee=zeros(P+1,P+1,nel);
%fee=zeros(P+1,nel);
Kzz=zeros(nel*P+1);
M=zeros(nel*P+1);
F=zeros(1,nel*P+1);
[x,w,H]=GetGLL(P+1);

for kk=1:nel
%c2=zeros(1,P+1);
weye=zeros(1,P+1);
%c2weye=zeros(1,P+1)
for ii=1:P+1
   %c2(ii)=vel(kk,x(ii),h)^2;
    weye(ii,ii)=w(ii);
    %c2weye(ii,ii)=w(ii)*c2(ii);
	%fee(ii,kk)=h(kk)/2*f(kk,ii,h)
end

%for ii=1:P+1
    %c2(ii)=vel(kk,x(ii),h)^2;
%end



mee(:,:,kk)=weye*h(kk)/2;

kzz(:,:,kk)=1/h(kk)*H*weye*H';

end




for kk=1:nel

    Kzz(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)=Kzz(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)+kzz(:,:,kk);
    M(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)=M(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)+mee(:,:,kk);
    %F(P*(kk-1)+1:P*(kk-1)+1+P)=F(P*(kk-1)+1:P*(kk-1)+1+P)+fee(:,kk);
end




d=zeros(1,nel*P+1);
dn1=zeros(1,nel*P+1);
v=zeros(1,nel*P+1);

dt=0.1002*0.85/2/2;
tn=40*10/dt;
t=0;
ntime=4;
c1=0.1786178958448091;
c2=-0.2123418310626054;
c3=-0.06626458266981849;

alpha=[ c1 c3 1-2*(c3+c1) c3 c1];
beta=[0.5-c2 c2 c2 0.5-c2];
f0=0.5;
t0=1.5/f0;


kkk=1;
tt=zeros(1,tn);
ft=zeros(1,tn);
dx0=zeros(1,tn);
for jj=1:tn
    t=t+dt;
dn1=d+dt*v;
 tau=(3.14*f0*(t-t0))^2; 
 F(nel*P/2+1)=(2*tau-1)*exp(-tau);
v=v+(dt*inv(M)*(-Kzz*dn1'+F'))';
d=dn1;

% for ii=1:ntime
% t=t+alpha(ii)*dt;
% d=d+alpha(ii)*dt*v;
% tau=(3.14*f0*(t-t0))^2;
% 
% F(nel*P/2+1)=(2*tau-1)*exp(-tau);
% v=v+(beta(ii)*dt*inv(M)*(-Kzz*d'+F'))';
% end
% 
% t=t+alpha(ntime+1);
% d=d+alpha(ntime+1)*v;
dx0(jj)=d(1);

tt(jj)=t;
ft(jj)=F(nel*P/2+1);


xx=linspace(0,1,nel*P+1);




    if mod(jj,1000)==0
        
        %subplot(8,4,kkk);
        figure(kkk);
        kkk=kkk+1;
        plot(xx,d,'r-*');
    end
end







figure(kkk);
plot(tt,dx0,'b-');
figure(kkk+1);
plot(tt,ft,'b-');















