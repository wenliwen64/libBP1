clear;
nel=100;
h=1/nel;
f=1/nel;
q0=1
T1=1;
fe=1:nel;
fe=fe*0+f;
fe(1)=q0+f/2;
fe(nel)=T1/h+f;
K=zeros(nel);
rou=1;
c=1;
M=eye(nel)*rou*c;
tn=150.000;
dt=0.005;
dn=fe*0;
ddn=fe*0;
deltaddn1=fe*0;
ntime=2;
alpha=0.5;
%





    
for ii=1:nel
    for jj=1:nel
        if ii==jj
            K(ii,jj)=2;
        end
        if ii==jj+1
            K(ii,jj)=-1;
        end
        if ii==jj-1
            K(ii,jj)=-1;
        end
    end
end
K(1,1)=1;
%K(100,100)=1;
K=K/h;
%d=inv(K)*fe';
for t=0:dt:tn
dn1=dn+dt*(1-alpha)*ddn;
ddn1=fe*0;
for ii=1:ntime
   R=fe-ddn1*M-dn1*K;
   deltaddn1=R*inv(M);
   dn1=dn1+dt*alpha*deltaddn1;
   ddn1=ddn1+deltaddn1;
end
dn=dn1;
ddn=ddn1;

figure(1);
xx=linspace(0,1,nel+1);




    if mod(t/dt,1000)==0
        
        subplot(8,4,t/dt/1000+1);
        plot(xx,[dn';1],'r');
    end
end
subplot(8,4,32);
x=linspace(0,1,100);
plot(x,T1+(1-x)*q0);
%plot(x,T1+(1-x)*q0+(1-x.^2)*0.5);
%title('f=0')
%legend('FEM(nel=10)','exact solution');
