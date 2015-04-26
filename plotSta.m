function plotSta(ret)
% function plotSta(ret)
% plot the stations with station name and political maps

if isfield(ret,'stacolor')
    stacolor=ret.stacolor;
else
    stacolor='b.';
end

load haiticoast.dat;
load west_hemi.dat;
load west_political.dat;
load worldcoast.dat
load ptimes;
load japan_coastline.dat;
station=['FUNV';'CRUV';'BAUV';'VIRV';'DABV';'MONV';'SIQV';'SANV';'QARV';'CURV';'CAPV';'SOCV';'CUPV';'JACV';'RIOV';'TURV';'TERV';'ELOV';'LUEV';'PCRV';'PRGV';'MAPV'];
 lat=[10.4692   10.6163    8.9433   10.5028   10.9218   11.9550   10.6488    9.5008   10.2065   10.0130    7.8647    8.2842   10.0563   11.0872    8.0618   10.4495    9.9637    7.0010 5.8432   10.1633    8.7600    9.8308];
 lon=-[66.8102   63.1842   68.0415   72.4060   70.6362   69.9703   69.8078   69.5363   70.5237   69.9612   72.3143   70.8563   65.7877   68.8298   61.8170   67.8395   69.1917   69.4833 61.4612   64.5897   64.6455   68.457];



%load /home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk;

[nel samlength]=size(ret.x);

r=ret.r;
lat0=ret.opr.lat0;
lon0=ret.opr.lon0;
stanm=ret.nm;


hold on;
% 
% plot(lon,lat,'r*');
plot(lon0,lat0,'r*');
plot( r(:,1), r(:,2),stacolor,'MarkerSize',20);
plot(west_hemi(:,1),west_hemi(:,2),'black');
plot(west_political(:,1),west_political(:,2),'black');
plot(japan_coastline(:,1),japan_coastline(:,2),'black');
plot(worldcoast(:,1),worldcoast(:,2),'black');
% xlim([-76 -64])
% ylim([6 16])
axis equal

for i=1:nel
%     text(r(i,1), r(i,2),[ stanm(i,:) ]);
%     text(r(i,1),r(i,2),ret.nm(i,:));
end
