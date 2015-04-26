clear;
close all;
nel=20;
h=[ones(1,nel/2)*75/50 ones(1,nel/2)*25/50]/100;
f=zeros(1,nel+1);
f(nel+1)=0;%??

diff=[ones(1,nel/2)*1 ones(1,nel/2)*1];
q0=1
T1=1;


rou=1;
c=1;


dt=0.00001;
tn=dt*300*40;
dn=zeros(1,nel);
ddn=zeros(1,nel);
deltaddn1=zeros(1,nel);
ntime=2;
alpha=0.5;
%
LM=zeros(2,nel);
for ii=1:2
    for jj=1:nel
        LM(ii,jj)=jj+ii-1;
    end
end
LM(2,nel)=0;
fee=zeros(2,nel);
kee=zeros(2,2,nel);
mee=zeros(2,2,nel);
K=zeros(nel);
M=zeros(nel);
fe=zeros(1,nel);
for kk=1:nel
fee(:,kk)=[f(kk)+2*f(kk+1) 2*f(kk)+f(kk+1)]'*h(kk)/6;
mee(:,:,kk)=eye(2)*rou*c*h(kk)/2;
kee(:,:,kk)=(eye(2)*2.-1)/h(kk)*diff(kk);
if kk==1
    fee(:,kk)=fee(:,kk)+[q0 0]';
end
if kk==nel
    fee(:,kk)=fee(:,kk)-[kee(1,2,kk) kee(2,2,kk)]'*T1;
end
end
for kk=1:nel-1
    K(kk:kk+1,kk:kk+1)=K(kk:kk+1,kk:kk+1)+kee(:,:,kk);
    M(kk:kk+1,kk:kk+1)=M(kk:kk+1,kk:kk+1)+mee(:,:,kk);
    fe(kk:kk+1)=fe(kk:kk+1)+fee(:,kk)';
end
K(nel,nel)=K(nel,nel)+kee(1,1,nel);
M(nel,nel)=M(nel,nel)+mee(1,1,nel);
fe(nel)=fe(nel)+fee(1,nel);
%M(1,1)=0.01;
%fe(nel)=100.001
%M=M*100;
kkk=1;
xx=zeros(1,nel+1);
xx(1)=0;
for iii=1:nel
    xx(iii+1)=h(iii)+xx(iii);
end
for t=0:dt:tn
dn1=dn+dt*(1-alpha)*ddn;
ddn1=zeros(1,nel);
for ii=1:ntime
   R=fe-ddn1*M-dn1*K;
   deltaddn1=R*inv(M);
   dn1=dn1+dt*alpha*deltaddn1;
   ddn1=ddn1+deltaddn1;
end
dn=dn1;
ddn=ddn1;

figure(1);





    if mod(t/dt,300)==0&&kkk<=31
        
        subplot(8,4,kkk);
        kkk=kkk+1;
        plot(xx,[dn 1],'r-*');
        xlim([0 xx(nel+1)]);
    end
end
subplot(8,4,32);
x=linspace(0,1,100);
plot(x,T1+(1-x)*q0);

%plot(x,T1+(1-x)*q0+(1-x.^2)*0.5);
%title('f=0')
%legend('FEM(nel=10)','exact solution');