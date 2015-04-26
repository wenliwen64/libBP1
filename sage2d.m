clear all;
M=2;%number of signals
k=@(lamda,theta) 2*pi/lamda*[cos(theta) sin(theta)];% wave vector
v=@(r,k) exp(1i*(r*k'));%steering vector
delay= @(x,y,the,c) (x/cos(the)+(y-tan(the)*x)*sin(the))/c;
d=0.05;%km
nel=14;
%km
r=zeros(nel,2);%coordinates of the array element
%r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
r(:,1)=[-371 -298 -192 38 -23 0 23 177 213 260 98 -72 -131 -195]/1000;
r(:,2)=[-302 -208 -284 -180 -13 0 -9 98 235 411 213 328 393 356]/1000;
time=25;%total time in seconds
samplerate=200;
t=linspace(0,time,time*samplerate);
L=time*samplerate;
fc=30;
f=zeros(1,M);
f(1)=2;%Hz
f(2)=3.0;
f(3)=4.0;
vp=zeros(1,M);
vp(1)=5;%km
vp(2)=4;
vp(3)=3;
phic=zeros(1,M);
phic(1)=pi/4;
phic(2)=pi/2;
phic(3)=pi;
SNR=10000;
theta=zeros(1,M);
theta(1)=120*pi/180;
theta(2)=110*pi/180;
theta(3)=130*pi/180;

% S=zeros(1,samplerate*time);
% for i=1:samplerate *time
%     if t>=0
% S(i)=A*(sin(2*pi*f0*t(i))+sin(2*pi*1.2*t(i)-2)+sin(2*pi*0.8*t(i)+2));
%     end
% end
S=zeros(M,L);
X=zeros(nel,L);
U=zeros(nel,L);
for i=1:nel
     U(i,:)=hilbert(1/SNR*randn(1,time*samplerate));
end
H=zeros(nel,M);
% for i=1:M
% S(i,:)=hilbert(sin(2*pi*t*f(i)).*sin(2*pi*fc*t+phic(i)));
% end
load y1;
load y2;
 S(1,:)=y2;
  S(2,:)=y1;
%  S(3,:)=y1;



H=zeros(nel,M);
for i=1:M
    H(:,i)=v(r,k(vp(i)/fc,theta(i)));
end
X=H*S+U;

Iteration=100;
s=zeros(M,L);
Z=zeros(M,L);
Rxx=zeros(nel,nel);
the=zeros(1,M);
the(1)=105*pi/180;
the(2)=101*pi/180;
the(3)=110*pi/180;
vpp=zeros(1,M);
vpp(1)=5.0;
vpp(2)=4.6;
vpp(3)=2.8;
cr=10;
cm=21;
vr=2;
vm=21;
for i=1:Iteration
    for m=1:M
        for q=1:M
            H(:,q)=v(r,k(vpp(q)/fc,the(q)));
        end
        
        Z=v(r,k(vpp(m)/fc,the(m)))*s(m,:)+X-H*s;
        Rxx=zeros(nel,nel);
        for p=1:L
            Rxx=Rxx+Z(:,p)*Z(:,p)';
        end
        R=Rxx/L;
        them=linspace(the(m)-cr*pi/180,the(m)+cr*pi/180,cm);
        vppm=linspace(vpp(m)-vr,vpp(m)+vr,vm);
        E=zeros(1,cm);
        for j=1:cm
           
            E(j)=v(r,k(vpp(m)/fc,them(j)))'*R*v(r,k(vpp(m)/fc,them(j)));
            
           
        end
        Emax=max(E);
        for j=1:cm
            
                if(E(j)==Emax)
                    %the(m)=them(j);
                    thenew=them(j);
                end
           
        end
        clear E;
        E=zeros(1,vm);
        for j=1:vm
           
            E(j)=v(r,k(vppm(j)/fc,the(m)))'*R*v(r,k(vppm(j)/fc,the(m)));
            
           
        end
        for j=1:vm
            
                if(E(j)==Emax)
                    vpp(m)=vppm(j);
                    
                end
           
        end
        the(m)=thenew;
        s(m,:)=1/nel*v(r,k(vpp(m)/fc,the(m)))'*Z;
        vn=1/nel*trace((eye(nel)-1/nel*v(r,k(vpp(m)/fc,the(m)))*v(r,k(vpp(m)/fc,the(m)))')*R);
        clear Z R E Emax;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    
    
    
    
    
    
    
end
    figure(1);
    for i=1:M
    subplot(M,1,i)            
    plot(t,S(i,:),'r',t,s(i,:),'b');
    xlim([1 3]);
    end
    
    
    
