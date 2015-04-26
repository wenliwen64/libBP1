clear;
close all;
load bulkj13;
bulk = bulkj13;
L=366*15-1550*2;
%L=366*15;

for k=1:10000
    sigma22(k)=(-bulk(k,2)+bulk(k,3))/2/(366*15-bulk(k,7)*2);
    sigma11(k)=(-bulk(k,4)+bulk(k,5))/2/2740;
    epsilon11(k)=(bulk(k,7)-1550)*2/L;
    %epsilon11(k)=(bulk(k,7))*2/L;
end
for k=0:499
    sig11(k+1)=sigma11(k*20+1);
    sig22(k+1)=sigma22(k*20+1);
    eps11(k+1)=epsilon11(k*20+1);
end
for k=1:499
    ds22(k)=sig22(k+1)-sig22(k);
    ds11(k)=sig11(k+1)-sig11(k);
    de11(k)=eps11(k+1)-eps11(k);
end
for k=1:499
    K(k)=(ds11(k)/6+ds22(k)/3)/de11(k);
    G(k)=(ds11(k)-ds22(k))/de11(k)/2;
end
for k=1:499
    K(k)=(sig11(k)/6+sig22(k)/3)/eps11(k);
    G(k)=(sig11(k)-sig22(k))/eps11(k)/2;
end
% for k=2:length(bulk)
%     K(k-1)=-(p(k)-p(k-1))/(v(k)-v(k-1))*v(k);
% end
figure(1);
hold on;
%plot(epsilon11,sigma22,'.');
%plot(epsilon11,sigma11,'r.');
plot(eps11,sig11,'r.',eps11,sig22,'b.');
xlim([0 inf]);
L=366*15-1550*2;
xlabel('strain');
ylabel('Sigma11(red) & Sigma22(blue) Pa');
figure(2);
hold on;
plot(eps11(2:500),K,'r');
plot(eps11(2:500),G,'b');
 
    xlim([0.0 inf]);
xlabel('strain');
ylabel('G(red) & K(blue) Pa');
%  ylim([-1e10 1e10])
% x=0:100:2e4;
% % a=1.5e6;
% b=-8e9;
% y=a*x+b;
% plot(x,y,'g')  ;
% kkk=0
% for kk=350000:50000:18800000
%     k=kk/1000+1;
%     plot(k,p(k),'r.');
%     kkk=kkk+1;
%     pp(kkk)=p(k);
%     vv(kkk)=v(k);
% end
% for k=2:length(pp)
%     K(k-1)=-(pp(k)-pp(k-1))/(vv(k)-vv(k-1))*v(k);
% end
% axis xy;
% figure(2);
% plot(vv(2:length(pp)),K,'-r.');
% %plot(bulk(2:length(bulk),1),K,'b-');
% %xlim([])
% 
% % for kk=30000:10000:160000
% %     k=kk/20+1;
% %     k1=(kk-10000)/20+1;
% %     K(kk/10000-2)=-(p(k)-p(k1))/(v(k)-v(k1));
% %     Kv(kk/10000-2)=vt(k);
% % end
% % figure;
% % plot(Kv,K,'*-')
% ylim([0 3e8])    
  