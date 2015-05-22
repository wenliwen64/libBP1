% this script will copy and shift the seismic timeseries for different stations
% according to a virtual epicenter. 
clear all;
load timeshift_newmoho.mat
load ./nsac4.mat
lat0 = ret.lat0;
lon0 = ret.lon0;
bparea_span = [-1., 1.; -1., 1.];

p0 = 11;
q0 = 31;
new_ep_lat = lat0 + bparea_span(1,1) + (bparea_span(1,2) - bparea_span(1,1))/40*p0
new_ep_lon = lon0 + bparea_span(2,1) + (bparea_span(2,2) - bparea_span(2,1))/40*q0

for i = 1:numel(ret.lat)
    ret.rdis(i)
    ret.sr*ret1.timeshift(q0, p0, i)
    ret.xori_new(i,:) = specshift(ret.xori(1,:), -ret.sr*ret1.timeshift(q0, p0, i));
%     plot(ret.xori_new(i,:), '*');
%     hold on;
%     plot(ret.xori(1,:), 'r');
end
ret.xori = ret.xori_new;
save('eu_testshift1.mat', 'ret', '-v7.3');
