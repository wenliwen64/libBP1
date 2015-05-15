%retracing test ...

%load /home/lwen/Documents/seis_research/raw_data/Hinet_nepal.v2.mat

ep_lat = ret.lat0;
ep_lon = ret.lon0;
grid_dim = [40 40];
bparea_span = [-5, 5; -5, 5];
lon_step = (bparea_span(1,2)-bparea_span(1,1)) / grid_dim(2);
lat_step = (bparea_span(2,2)-bparea_span(2,1)) / grid_dim(1);
evt_depth = 20;
loc_new_grid = zeros(grid_dim(1), grid_dim(2));
count = 1;
ii = 100;

new_grid = zeros(40, 40)
for i = 1:40 % latitude
    for j = 1:40 % longitude
        loc_grid = [-2+i*lat_step+ep_lat, -2+j*lon_step+ep_lon];
        loc_sta = [ret.lat(ii), ret.lon(ii)];
        loc_new_grid = retracing_liwen(loc_sta, loc_grid, evt_depth);
        loc_new_grid_y(count) = loc_new_grid(1);
        loc_new_grid_x(count) = loc_new_grid(2);
        loc_old_grid_x(count) = -2+j*lon_step+ep_lon;
        loc_old_grid_y(count) = -2+i*lat_step+ep_lat;
        
        new_grid_x(i, j) = loc_new_grid_x(count);
        new_grid_y(i, j) = loc_new_grid_y(count);
        new_grid_timeshift(i, j) = timeshift;        
        count = count+1;

        
    end
end

figure;
hold all;
scatter(loc_new_grid_x, loc_new_grid_y);
scatter(loc_old_grid_x, loc_old_grid_y);
