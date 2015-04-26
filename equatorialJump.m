figure(3);
 set(gcf,'Position',[10 10 600 700]);
subplot(311);
aaa=load('/export/scratch1/lsmeng/Ge161/matlab/orin10Nob10N140R50SP.dat'); 
a=aaa(:,5);
a=a-mean(a);
L=1;
H=1200;
a=a(L:H);
dis=aaa(L:H,1)+1400;
% a=hilshift(a,-5);
skew0=kurt(a,180,1);
a1=hilshift(a,skew0);
plot(dis,a);
title(['Skewness ' num2str(skew0) '^o']);
% xlabel('Distance (km)');
ylabel('Magnetic anomaly (nT)');
xlim([0 1000]);
ylim([-40 40]);

subplot(312);
aaa=load('/export/scratch1/lsmeng/Ge161/matlab/orin10Sob10N140R50SP.dat'); 
a=aaa(:,5);
a=a-mean(a);
L=1;
H=1200;
a=a(L:H);
dis=aaa(L:H,1)+1400;
% a=hilshift(a,-5);
skew0=kurt(a,180,1);
a1=hilshift(a,skew0);
plot(dis,a);
title(['Skewness ' num2str(skew0) '^o']);
% xlabel('Distance (km)');
ylabel('Magnetic anomaly (nT)');
xlim([0 1000]);
ylim([-40 40]);

subplot(313);
plot(dis,a1);
title(['Skewness ' num2str(0) '^o']);
xlim([0 1000]);
ylim([-40 40]);
xlabel('Distance (km)');
ylabel('Magnetic anomaly (nT)');
