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
d=inv(K/h)*fe';
xx=linspace(0,1,nel+1)
plot(xx,[d;1],'r*');
hold on;
x=linspace(0,1,100);
%plot(x,T1+(1-x)*q0);
plot(x,T1+(1-x)*q0+(1-x.^2)*0.5)
title('f=0')
legend('FEM(nel=10)','exact solution');