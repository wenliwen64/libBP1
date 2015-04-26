function dep=intplaneDip(lon,lat)

uz=utmzone( lat,lon);
uz=str2num(uz(1:2));
[X,Y,zone]=wgs2utm(lat,lon,uz,'N'); 
X=X/1e3;
Y=Y/1e3;

Y0=6.0949e+03;X0=531.9841;%(55N, 153.5E)
Y1=5.9835e+03;X1=486.8901;%(54N, 152.8E)
% Y1=Y1-Y0;X1=X1-X0;
% Y0=0;X0=0;

a=1/(X1-X0);
b=-1/(Y1-Y0);
c=-X0/(X1-X0)+Y0/(Y1-Y0);
d=sqrt(a^2+b^2);

% a = -0.0090;
% b = 0.0222;
% c = 42.9146;
% d = 0.0239;

% theta=atan(30/90);


dep=600-(a*X+b*Y+c)/d;

% dep=620+D*30/90;
% Z0=43.3631;
% dep=-interp2(x,y,z,X,Y)-Z0;