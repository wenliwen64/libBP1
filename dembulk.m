clear;
close all;
load bulk20;
bulk = bulk20;
for k=1:length(bulk)
    p(k)=(-bulk(k,2)+bulk(k,3)-bulk(k,4)+bulk(k,5))/4/(880-2*bulk(k,8))^2;
    v(k)=(880-2*bulk(k,8))^3;
    vt(k)=1-v(k)/860^3;
end
% for k=2:length(bulk)
%     K(k-1)=-(p(k)-p(k-1))/(v(k)-v(k-1))*v(k);
% end
figure(1);
hold on;
plot(p,'.')
kkk=0
for kk=350000:50000:19950000
    k=kk/1000+1;
    plot(k,p(k),'r.');
    kkk=kkk+1;
    pp(kkk)=p(k);
    vv(kkk)=v(k);
end
for k=2:length(pp)
    K(k-1)=-(pp(k)-pp(k-1))/(vv(k)-vv(k-1))*v(k);
end
axis xy;
figure(2);
plot(vv(2:length(pp)),K,'-r.');
%plot(bulk(2:length(bulk),1),K,'b-');
%xlim([])

% for kk=30000:10000:160000
%     k=kk/20+1;
%     k1=(kk-10000)/20+1;
%     K(kk/10000-2)=-(p(k)-p(k1))/(v(k)-v(k1));
%     Kv(kk/10000-2)=vt(k);
% end
% figure;
% plot(Kv,K,'*-')
ylim([0 3e8])    
    