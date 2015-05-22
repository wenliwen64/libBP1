% this script will copy and shift the seismic timeseries for different stations
% according to a virtual epicenter. 
clear all;
load timeshift_newmoho.mat
load ./nsac4.mat
lat0 = ret.lat0;
lon0 = ret.lon0;

p0 = 10;
q0 = 30;
new_ep_lat = lat0 - 1.4 + 2.8/40*p0;
new_ep_lon = lon0 - 1.4 + 2.8/40*q0;

for i = 1:numel(ret.lat)
    ret.rdis(i)
    ret.sr*ret1.timeshift( p0, q0,i)
    ret.xori_new(i,:) = specshift(ret.xori(1,:), ret.sr*ret1.timeshift(p0, q0, i));
end
ret.xori = ret.xori_new;
save('eu_testshift1.mat', 'ret', '-v7.3');