% retracing test ...
% You have to load a ret data first, and then run this script.
% load /home/lwen/Documents/seis_research/raw_data/Hinet_nepal.v2.mat
javaaddpath /home/lwen/Documents/seis_research/libBP1/taup/matTaup.jar
ep_lat = ret.lat0;
ep_lon = ret.lon0;
grid_dim = [40 40];
bparea_span = [-1.4, 1.4; -1.4, 1.4];
lon_step = (bparea_span(1,2)-bparea_span(1,1)) / grid_dim(2);
lat_step = (bparea_span(2,2)-bparea_span(2,1)) / grid_dim(1);
evt_depth = 20;
loc_new_grid = zeros(grid_dim(1), grid_dim(2));
count = 1;
ii = 1;
load ptimes.mat
new_grid = zeros(40, 40);
nsta = numel(ret.lat);
new_grid_timeshift = zeros(nsta, 40, 40);
old_grid_timeshift = zeros(nsta, 40, 40);
for ii = 1:nsta
    ii
    [dummy, timeshift0] = retracing_liwen([ret.lat(ii), ret.lon(ii)], [ep_lat, ep_lon], evt_depth);
    count = 1;
    for i = 1:40 % latitude
        for j = 1:40 % longitude
            loc_grid = [-1.4+i*lat_step+ep_lat, -1.4+j*lon_step+ep_lon];
            loc_sta = [ret.lat(ii), ret.lon(ii)];
            [loc_new_grid, timeshift_temp] = retracing_liwen(loc_sta, loc_grid, evt_depth);
            loc_new_grid_y(count) = loc_new_grid(1);
            loc_new_grid_x(count) = loc_new_grid(2);
            loc_old_grid_x(count) = -2+j*lon_step+ep_lon;
            loc_old_grid_y(count) = -2+i*lat_step+ep_lat;
            sd = phtime(ret.fl, ep_lat, ep_lon, loc_grid(1), loc_grid(2), loc_sta(1), loc_sta(2), rr, tt);
            timeshift(count) = timeshift_temp+sd-timeshift0;
            timeshift_temp;% + sd;
            timeshift_old(count) = sd;

            new_grid_x(i, j) = loc_new_grid_x(count);
            new_grid_y(i, j) = loc_new_grid_y(count);
            new_grid_timeshift(ii, i, j) = timeshift(count);
            old_grid_timeshift(ii, i, j) = timeshift_old(count);
            

            count = count + 1;
        end
    end
    [gridx, gridy] = ndgrid(ep_lon+bparea_span(1,1):lon_step:ep_lon+bparea_span(1,2), ep_lat+bparea_span(2,1):lat_step:ep_lat+bparea_span(2,2));
    F1 = scatteredInterpolant(loc_new_grid_x(:), loc_new_grid_y(:), timeshift(:));
    F2 = scatteredInterpolant(loc_old_grid_x(:), loc_old_grid_y(:), timeshift_old(:));

    v1 = F1(gridx, gridy);
    v2 = F2(gridx, gridy);
    
    ret1.timeshift(ii, :, :) = new_grid_timeshift(ii, :, :);
    ret1.stalat(ii) = ret.lat(ii);
    ret1.stalon(ii) = ret.lon(ii);
end
ret.lat_step = lat_step;
ret.lon_step = lon_step;

save('timeshift_newmoho.mat', 'ret1', '-v7.3');

% figure;
% hold all;
% scatter(loc_new_grid_x, loc_new_grid_y);
% scatter(loc_old_grid_x, loc_old_grid_y);
% 
% [gridx, gridy] = ndgrid(ep_lon+bparea_span(1,1):lon_step:ep_lon+bparea_span(1,2), ep_lat+bparea_span(2,1):lat_step:ep_lat+bparea_span(2,2));
% F1 = scatteredInterpolant(loc_new_grid_x(:), loc_new_grid_y(:), timeshift(:));
% F2 = scatteredInterpolant(loc_old_grid_x(:), loc_old_grid_y(:), timeshift_old(:));
% 
% v1 = F1(gridx, gridy);
% v2 = F2(gridx, gridy);
% figure;
% 
% surf(gridx, gridy, v1);
% colorbar;
% figure;
% surf(gridx, gridy, v2);
% colorbar;

