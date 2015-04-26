function [Pm,Uv,A]=musicDOA(yy,time,tts1)
M=1;%number of signals
k=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.3^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
delay= @(x,y,the,c) (x/cos(the)+(y-tan(the)*x)*sin(the))/c;
d=0.05;%km
nel=12;
%km
r=zeros(nel,3);%coordinates of the array element
%r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
% r(:,1)=[-371 -298 -192 38 -23 0 23 177 213 260 98 -72 -131 -195]/1000;
% r(:,2)=[-302 -208 -284 -180 -13 0 -9 98 235 411 213 328 393 356]/1000;
% r(:,3)=[-25 -24 -16 -18 -4 0 2 18 12 3 -1 -16 -4 -2]/1000;
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
%time=2;%total time in seconds
samplerate=200;
t=linspace(0,time,time*samplerate);
L=time*samplerate;
fc=11;
f=zeros(1,M);
f(1)=2.2;%Hz
f(2)=1.0;
f(3)=3;
vp=0.9;

phic=zeros(1,M);
phic(1)=pi/4;
phic(2)=pi/2;
phic(3)=pi;
SNR=1;
theta=zeros(1,M);
phi=zeros(1,M);
theta(1)=120*pi/180;
theta(2)=110*pi/180;
theta(3)=130*pi/180;
phi(1)=-10*pi/180;
phi(2)=-20*pi/180;
phi(3)=-30*pi/180;
ux=[0.2 0.4 0.6];
uy=[0.3 0.4 0.1];

% S=zeros(1,samplerate*time);
% for i=1:samplerate *time
%     if t>=0
% S(i)=(sin(2*pi*fc*t(i))+sin(2*pi*1.2*t(i)-2)+sin(2*pi*0.8*t(i)+2));
%     end
% end
S=zeros(M,L);
X=zeros(nel,L);
U=zeros(nel,L);
for i=1:nel
     U(i,:)=hilbert(1/SNR*std(yy(1,:))*randn(1,L));
end
H=zeros(nel,M);
for i=1:M
%S(i,:)=hilbert(sin(2*pi*t*f(i)).*sin(2*pi*fc*t+phic(i)));
%S(i,:)=hilbert(sin(2*pi*(fc+f(i))*t+phic(i)));
 %S(i,:)=hilbert(sin(2*pi*t*f(i)).*yy(1,:));
 S(i,:)=hilbert(yy(1,:));
end
% load y1;
% load y2;
% load y3;
%  S(1,:)=y1;
%   S(2,:)=y2;
%   S(3,:)=y3;
%load x;

for i=1:nel
    X(i,:)=hilbert(yy(i,:));
end

% 
%  H=zeros(nel,M);
% for i=1:M
%     H(:,i)=v(r,k(fc,ux(i),uy(i)));
% end
% X=H*S+U;

Iteration=80;
s=zeros(M,L);
Z=zeros(M,L);
Rxx=zeros(nel,nel);
the=zeros(1,M);
the(1)=-100*pi/180;
the(2)=-30*pi/180;
the(3)=-30*pi/180;
ph=zeros(1,M);
ph(1)=-64*pi/180;
ph(2)=-10*pi/180;
ph(3)=-10*pi/180;
cr=15;
cm=21;
vr=15;
vm=21;
%Pm=zeros(360,181);
for p=1:L
    Rxx=Rxx+X(:,p)*X(:,p)';
end
[Uv,A]=eig(Rxx);
un=zeros(nel,nel);
un(:,1:nel-5)=Uv(:,1:nel-5)
Un=un*un'
ps=100;
qs=50;
X=zeros(ps,qs);
Y=zeros(ps,qs);

ux=linspace(-1,1,ps);
uy=linspace(-1,1,qs);
for q=1:qs
    X(:,q)=ux;
end
for p=1:ps
    Y(p,:)=uy';
end
for p=1:ps
    for q=1:qs
       
        clear a;
        a=v(r,k(fc,ux(p),uy(q)));
        Pm(p,q)=a'*a/(a'*Un*a);
    end
end

figure(10);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
contourf(X,Y,real(Pm),15);
%caxis([0 15]);
colorbar;
title(['window start at' num2str(tts1) 'sec']);
%shading interp
% for i=1:Iteration
%     for m=1:M
%         for q=1:M
%             H(:,q)=v(r,k(vp/fc,the(q),ph(q)));
%         end
%         
%         Z=v(r,k(vp/fc,the(m),ph(m)))*s(m,:)+X-H*s;
%         Rxx=zeros(nel,nel);
%         for p=1:L
%             Rxx=Rxx+Z(:,p)*Z(:,p)';
%         end
%         R=Rxx/L;
%         them=linspace(the(m)-cr*pi/180,the(m)+cr*pi/180,cm);
%         phm=linspace(ph(m)-vr*pi/180,ph(m)+vr*pi/180,vm);
%         E=zeros(cm,vm);
%         for j=1:cm
%              for n=1:vm
%             E(j,n)=v(r,k(vp/fc,them(j),phm(n)))'*R*v(r,k(vp/fc,them(j),phm(n)));
%             
%              end
%         end
%         Emax=max(max(E));
%         for j=1:cm
%             for n=1:vm
%                 if(E(j,n)==Emax)
%                     the(m)=them(j);
%                       ph(m)=phm(n);
%                 end
%             end
%         end
%       
%         
%         s(m,:)=1/nel*v(r,k(vp/fc,the(m),ph(m)))'*Z;
%         vn=1/nel*trace((eye(nel)-1/nel*v(r,k(vp/fc,the(m),ph(m)))*v(r,k(vp/fc,the(m),ph(m)))')*R);
%         clear Z R E Emax;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%     end
%     
%     
%     
%     
%     
%     
%     
%     
% end
% 
% for q=1:M
%     H(:,q)=v(r,k(vp/fc,the(q),ph(q)));
% end
%         Er=X-H*s;
%         MEr=0;
% for p=1:L
%         MEr=MEr+Er(:,p)'*Er(:,p);
% end
      
%     figure(1);
%     for i=1:M
%     subplot(M,1,i)            
%     plot(t,S(i,:),'r',t,s(i,:),'b');
%     xlim([0 time]);
%     end
%     
    
