function plotStatele(DataFile)
load(DataFile);
load worldcoast.dat
figure(40);
close(40);
figure(40);
hold on;
r=ret.r;
lat0=ret.lat0;
lon0=ret.lon0;
hold on;
plot(lon0,lat0,'r*');
plot( r(:,1), r(:,2),'b.','MarkerSize',20);
plot(worldcoast(:,1),worldcoast(:,2),'black');
axis equal
