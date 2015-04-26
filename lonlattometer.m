function [x y]=lonlattometer(lons,lats)
lonm=mean(lons);
latm=mean(lats);
lonm=lons(1);
latm=lats(1);
x=110.49*(lons-lonm);
y=110.49*cosd(latm)*(lats-latm);