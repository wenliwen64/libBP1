% % % %%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
load movieBP;
Power=Power/max(Power(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nn mm]=size(Power);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h12=figure(12);
set(h12, 'Position', [100, 100, 600, 600]);
set(gcf, 'PaperPositionMode', 'auto') ;

cmm=colormap;
hold on;
c0=linspace(1,64,tend);
r=interp1(1:64,cmm(:,1),c0);
g=interp1(1:64,cmm(:,2),c0);
b=interp1(1:64,cmm(:,3),c0);
r1=interp1(1:tend,r,t);
g1=interp1(1:tend,g,t);
b1=interp1(1:tend,b,t);
kt=0;
for j=1:length(t)
    kt=kt+1;
    for jj=1:1
        marker='o';
        scatter(buy(j,jj),bux(j,jj),abs(Power(j,jj))*150,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
axis equal;
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Se (s/km)');
ylabel('Sn (s/km)');
box on;
hold off;
print('-dpdf','-r300','summaryBP.pdf') 
x3=[ t' bux(:,1) buy(:,1) Power(:,1)];
save('Slownessdots','x3','-ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h16=figure(16);
imagesc(dist,ret.begin:ret.step:ret.end,rup);
axis xy;
xlabel('Distance (km)');
ylabel('Time (s)');
box on
colorbar;
% caxis([0 0.7]);
shading flat;
hold;
for i=1:length(maxcc)
    scatter(distpeak(i),t(i),max([maxcc(i)-cct 0.0001])*200,'white','o','filled');
    scatter(distpeak(i),t(i),max([maxcc(i)-cct 0.0001])*200,'red','o');
end
 distpeakL=distpeak(maxcc>cct);
 tL=t(maxcc>cct);
 maxccL=maxcc(maxcc>cct);
 x4=[ tL' distpeakL' maxccL'];
save('Distancedots','x4','-ascii');
chsize(15);
set(h16, 'Position', [100, 100, 600, 400]); 
set(gcf, 'PaperPositionMode', 'auto') ;
print('-depsc','-r600','powervstime.eps') 
print('-dpdf','-r600','powervstime.pdf') 
print('-djpeg','-r600','powervstime.jpg') 
  