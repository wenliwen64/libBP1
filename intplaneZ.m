function dep=intplaneZ(lon,lat)
load fplane;
uz=utmzone( lat,lon);
uz=str2num(uz(1:2));
[X,Y,zone]=wgs2utm(lat,lon,uz,'N'); 
X=X/1e3;
Y=Y/1e3;
Z0=43.3631;
dep=-interp2(x,y,z,X,Y)-Z0;