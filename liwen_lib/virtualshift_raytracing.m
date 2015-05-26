% this script will copy and shift the seismic timeseries for different stations
% according to a virtual epicenter. 
clear all;
close all;
load timeshift_newmoho_eu.mat
load ./nsac4.mat
lat0 = ret.lat0;
lon0 = ret.lon0;
bparea_span = [-1., 1.; -1., 1.];
xslices = 40;
yslices = 40;
p0 = 15; % y
q0 = 25; % x
bpgrid_x_vec = linspace(lon0 + bparea_span(1, 1), lon0 + bparea_span(1,2), xslices);
bpgrid_y_vec = linspace(lat0 + bparea_span(2, 1), lat0 + bparea_span(2, 2), yslices);
%new_ep_lat = bpgrid_x_vec(q0);lat0 + bparea_span(1,1) + (bparea_span(1,2) - bparea_span(1,1))/40*p0
%new_ep_lon = lon0 + bparea_span(2,1) + (bparea_span(2,2) - bparea_span(2,1))/40*q0
figure
ret.x = ret.xori;
plotAll1(ret);
for i = 1:numel(ret.lat)
    ret.rdis(i);
    ret.sr*ret1.timeshift(q0, p0, i)
     ret.xori_new(i,:) = specshift(ret.xori(2,:), -ret.sr*ret1.timeshift(q0, p0, i));
    %ret.xori_new(i,:) = specshift(ret.xori(2,:), 0);
%     plot(ret.xori_new(i,:), '*');
%     hold on;
%     plot(ret.xori(1,:), 'r');
end
ret.xori = ret.xori_new;
ret.x = ret.xori_new;
figure;
plotAll1(ret);
save('eu_testshift1.mat', 'ret', '-v7.3');
