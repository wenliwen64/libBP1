function [xr xt]=rotateseis(lat0,lon0,lat,lon,xn,xe)

phi=azimuth(lat0,lon0,lat,lon)-90;
xr=sind(phi)*xn+cosd(phi)*xe;
xt=cosd(phi)*xn+sind(phi)*xe;