% % % %%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
load movieBP;
load kurslabcontourmat;
% load ../earlyaftershock;
tend=50;
Power=Power/max(Power(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    for jj=1:min([ nw(kt) 1])
        marker='o';
        scatter(buy(kt,jj),bux(kt,jj),abs(Power(kt,jj))*150,[r1(j) g1(j) b1(j)],marker,'filled');
    end
end
axis equal;
plot(worldcoast(:,1),worldcoast(:,2),'black','LineWidth',2);
plot(lon0,lat0,'r*','MarkerSize',10);
plot(kurslabcontourmat(:,1),kurslabcontourmat(:,2),'black','LineWidth',2);
% plot(earlyaftershock(:,2),earlyaftershock(:,1),'k.');

latC=-19.85;
lonC=-70.7;

plot(lonC,latC,'b*','MarkerSize',10);

set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1])
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);
% xlim([-71.2 -70.4]);
% ylim([-20.4 -19]);
colorbar;
colormap(cmm);
caxis([0 tend]);
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
box on;
x3=[ t' bux(1:length(t),1) buy(1:length(t),1) Power(1:length(t),1)/max(Power(1:length(t),1))];
save('HFdots','x3','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0 = 100;
uz0=utmzone(lat0,lon0);
uz=str2num(uz0(1:2));
if lat0<0
NS = 'S';
else
NS = 'N';
end
[X0,Y0,zone]=wgs2utm(lat0,lon0,uz,NS);
X0=X0/1e3;
Y0=Y0/1e3;
X1=X0+100*cosd(theta0);
Y1=Y0+100*sind(theta0);
[lat1 lon1]=utm2wgs(X1*1e3,Y1*1e3,uz0);
% plot([lon0 lon1],[lat0 lat1],'r');
title('0.5 - 2 s')
hold off;
print('-dpdf','-r300','summaryBP.pdf') 
%%%%%%%%%%%%%%%%%%%%%%%
mp1=max(Power(:,1));
mp2=max(Power(:,2));

R0 = [ [cosd(theta0) sind(theta0)] ; [-sind(theta0) cosd(theta0)] ];


kt=0;
xymax=100;
discmm=linspace(-xymax,xymax,length(cmm(:,1)));



h20=figure(20);
  set(h20, 'Position', [100, 100, 500, 500]); 
  set(gcf, 'PaperPositionMode', 'auto') ;

  coor0 = [lon0 lat0];
   coor10 = (R0*coor0')';

  
for j=1:step:tend
    kt=kt+1;

    for jj=1:min([ nw(kt) 2])
        
        
        coor = [buy(kt,jj) bux(kt,jj)];
            coor1 = (R0*coor')';
            xx=(coor1(1)-coor10(1))*110.49;
            yy=(coor1(2)-coor10(2))*110.49;

        
        rx=interp1(discmm,cmm(:,1),xx);
        gx=interp1(discmm,cmm(:,2),xx);
        bx=interp1(discmm,cmm(:,3),xx);
        
        ry=interp1(discmm,cmm(:,1),yy);
        gy=interp1(discmm,cmm(:,2),yy);
        by=interp1(discmm,cmm(:,3),yy);
        subplot(211);
        hold on;
        if jj==1
            scatter(xx,t(j),Power(kt,jj)*50/mp1,'g','o','filled');
        else
           % scatter(xx,t(j),Power(kt,jj)*20/mp2,'g','s');
        end
               subplot(212);
        hold on;
        if jj==1
            scatter(yy,t(j),Power(kt,jj)*50/mp1,'g','o','filled');
        else
           % scatter(xx,t(j),Power(kt,jj)*20/mp2,'g','s');
        end
       
    end

end
subplot(211);
xlabel('Distance along X (km)');
ylabel('Time (s)');
xlim([-xymax xymax]);
ylim([0 tend]);
chsize(15);
box on;
subplot(212);
xlabel('Distance along Y (km)');
ylabel('Time (s)');
xlim([-xymax xymax]);
ylim([0 tend]);
chsize(15);
box on;
 print('-dpdf','-r300','distXtYt.pdf') ;
 
 %%%%%%%%%%%%%%%%
 figure(21)
 plot(t,Power(1:length(t),1),'r-.');
 print('-dpdf','-r300','Power.pdf') ;
