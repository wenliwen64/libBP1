

clear all;
close all;
%case A 
nel=24;
P=5;
omiga=2*3.14*2.5;
h=[ones(1,8)*300/4 ones(1,16)*850/4];







kxx=zeros(P+1,P+1,nel);
kzz=zeros(P+1,P+1,nel);
mee=zeros(P+1,P+1,nel);
Kxx=zeros(nel*P+1);
Kzz=zeros(nel*P+1);
M=zeros(nel*P+1);

[x,w,H]=GetGLL(P+1);

for kk=1:nel
c2=zeros(1,P+1);
weye=zeros(1,P+1);
c2weye=zeros(1,P+1)
for ii=1:P+1
    c2(ii)=vel(kk,x(ii),h)^2;
    weye(ii,ii)=w(ii);
    c2weye(ii,ii)=w(ii)*c2(ii);
end

for ii=1:P+1
    c2(ii)=vel(kk,x(ii),h)^2;
end



mee(:,:,kk)=weye*h(kk)/2;
kxx(:,:,kk)=c2weye*h(kk)/2;
kzz(:,:,kk)=1/h(kk)*H*c2weye*H';

end




for kk=1:nel
    Kxx(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)=Kxx(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)+kxx(:,:,kk);
    Kzz(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)=Kzz(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)+kzz(:,:,kk);
    M(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)=M(P*(kk-1)+1:P*(kk-1)+1+P,P*(kk-1)+1:P*(kk-1)+1+P)+mee(:,:,kk);
end

%K(nel,nel)=K(nel,nel)+kee(1,1,nel);
%M(nel,nel)=M(nel,nel)+mee(1,1,nel);
%fe(nel)=fe(nel)+fee(1,nel);
%M(1,1)=0.01;
%fe(nel)=100.001
%M=M*100;
Kzz=0.5*(Kzz+Kzz');
TA=zeros(nel*P+1);
TB=zeros(nel*P+1);
TA=Kzz-omiga^2*M;
TB=Kxx;
A=TA(1:nel*P,1:nel*P);
B=TB(1:nel*P,1:nel*P);
[um km]=eigs(A,B,2,'SA');

jj=1;
for ii=1:nel
    for kk=1:P
        zm(jj)=dep(ii,x(kk),h);
        jj=jj+1;
    end
end
[u1 k]=love_analytic(zm,0);
[u2 k]=love_analytic(zm,1);
umn1=um(:,1)/um(1,1);
umn2=um(:,2)/um(1,2);
%plot(zz,u1,'b',zz,u2,'r',zm,(um(:,1)-min(um(:,1)))/abs(max(um(:,1))-min(um(:,1))),'b*',zm,(um(:,2)-min(um(:,2)))/abs(max(um(:,2))-min(um(:,2)))-1,'b*');
plot(zm,u1,'b',zm,u2,'r',zm,umn1,'bo',zm,umn2,'r*');
legend('analytic for mode 0','analytic for mode 1','SEM for mode 0','SEM for mode 1');
title('love wave : um vs depth')
error1=0;
error2=0;
for ii=2:size(um(:,1))
error1=error1+(zm(ii)-zm(ii-1))*(umn1(ii)-u1(ii))^2;
error2=error2+(zm(ii)-zm(ii-1))*(umn2(ii)-u2(ii))^2;
end
ylim([-1 1]);
error1
error2
km


figure(2);

ha=[6 12 24];
ea1=[1.8384 0.1598 0.3411];
ea2=[14.0484 1.8587 6.2527];

hb=[5 10 20];
eb1=[0.8693 0.3083 0.1646];
eb2=[10.0361 5.8903 4.5515];
semilogy(4000./ha,ea1,'r-o',4000./ha,ea2,'r-*',4000./hb,eb1,'b-o',4000./hb,eb2,'b-*');
legend('mode0 caseA','mode1 caseA','mode0 caseB','mode1 caseB');
title('error as a function of element size')

