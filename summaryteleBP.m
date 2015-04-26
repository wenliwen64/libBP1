% % % %%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
load movieBP;
load kurslabcontourmat;
Power=Power/max(Power(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tend=50;
t=1:step:tend;
dep0=20;
load ptimes;
[nn mm]=size(Power);
for j=1:length(t)
    sd=phtime(1,lat0,lon0,bux(j,1),buy(j,1),ret.r(:,2),ret.r(:,1),rr,tt)'; 
    t(j)=t(j)-mean(sd);
end

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
        scatter(buy(kt,jj),bux(kt,jj),abs(Power(kt,jj))*150,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
axis equal;
plot(worldcoast(:,1),worldcoast(:,2),'black','LineWidth',2);
plot(lon0,lat0,'r*','MarkerSize',10);
plot(kurslabcontourmat(:,1),kurslabcontourmat(:,2),'black','LineWidth',2);

set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1])
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
box on;
hold off;
print('-dpdf','-r300','summaryBP.pdf') 
x3=[ t' bux(:,1) buy(:,1) Power(:,1)/max(Power(:,1))];
save('HFdots','x3','-ascii');