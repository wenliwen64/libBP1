function dep=intplaneZC(lon,lat)
load('~/wk/matlab/bp/russia82/uzzcurve.mat');
dep=interp2(uy,ux,uz,lon,lat);