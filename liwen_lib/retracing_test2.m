% retracing test ...
% You have to load a ret data first, and then run this script.
% load /home/lwen/Documents/seis_research/raw_data/Hinet_nepal.v2.mat

load nsac4.mat
d2km = 6371*2*3.1415926/360;
km2d = 1/d2km;
%load nsac5_au.mat
ep_lat = ret.lat0;
ep_lon = ret.lon0;
epicenter = [ep_lat, ep_lon];
xslices = 40;
yslices = 40;
grid_dim = [xslices yslices];
bparea_span = [-1., 1.; -1., 1.];
lon_step = (bparea_span(1,2)-bparea_span(1,1)) / (grid_dim(2)-1);
lat_step = (bparea_span(2,2)-bparea_span(2,1)) / (grid_dim(1)-1);
evt_depth = 20;
loc_new_grid = zeros(grid_dim(1), grid_dim(2));
count = 1;
ii = 1;
load ptimes.mat
new_grid = zeros(xslices, yslices);
nsta = numel(ret.lat);
new_grid_timeshift = zeros(xslices, yslices, nsta);
old_grid_timeshift = zeros(xslices, yslices, nsta);
loc_new_grid_x = zeros(40*40, 1);
loc_flat_grid_x = zeros(40*40, 1);
loc_new_grid_y = zeros(40*40, 1);
loc_flat_grid_y = zeros(40*40, 1);
loc_old_grid_x = zeros(40*40, 1);
loc_old_grid_y = zeros(40*40, 1);
timeshift = zeros(40*40, 1);
timeshift_flat = zeros(40*40, 1);
new_grid_x = zeros(40, 40);
%ret1.timeshift = zeros(41, 41, nsta);
bpgrid_x_vec = linspace(ep_lon + bparea_span(1, 1), ep_lon + bparea_span(1,2), xslices);
bpgrid_y_vec = linspace(ep_lat + bparea_span(2, 1), ep_lat + bparea_span(2, 2), yslices);


for ii = 1:1%nsta
    ii
%     ep_retracing_input.loc_sta = [ret.lat(ii), ret.lon(ii)];
%     ep_retracing_input.loc_grid = [ep_lat, ep_lon];
%     ep_retracing_input.evt_depth = 20;
%     ep_retracing_input.third_layer_depth_retracing = 20;
%     ep_retracing_input.start_depth = 80;
%     ep_retracing_input.epicenter_r = [ep_lon, ep_lat]*d2km;
%     ep_retracing_input.bprange_r = bparea_span*d2km;
%     
%     ep_retracing_output = retracing_liwen2(ep_retracing_input);
%     
    count = 1;
    for i = 1:yslices % latitude
        for j = 1:xslices % longitud
            
            loc_grid = [bpgrid_y_vec(i), bpgrid_x_vec(j)];%[bparea_span(1,1)+(i-1)*lat_step+ep_lat, bparea_span(2,1)+(j-1)*lon_step+ep_lon];  % lat-lon system
            loc_sta = [ret.lat(ii), ret.lon(ii)];
            input_temp.loc_sta = loc_sta;% lat-lon system
            input_temp.loc_grid = loc_grid;
            input_temp.evt_depth = 20;
            input_temp.third_layer_depth_retracing = 20;
            input_temp.start_depth = 80;
            input_temp.epicenter_r = [ep_lon, ep_lat]*d2km;
            input_temp.bprange_r = bparea_span*d2km;
            
            output_temp = retracing_liwen2(input_temp);
            
            loc_new_grid_x(count) = output_temp.loc_new_grid(2);
            loc_new_grid_y(count) = output_temp.loc_new_grid(1);
            loc_flat_grid_x(count) = output_temp.loc_new_grid_flat(2);
            loc_flat_grid_y(count) = output_temp.loc_new_grid_flat(1);
            loc_old_grid_x(count) = loc_grid(2); % bparea_span(2,1) + (j-1)*lon_step + ep_lon;
            loc_old_grid_y(count) = loc_grid(1); % bparea_span(1,1) + (i-1)*lat_step + ep_lat;
            %sd = phtime(ret.fl, ep_lat, ep_lon, loc_new_grid(1), loc_new_grid(2), loc_sta(1), loc_sta(2), rr, tt);
            %sd = phtime(ret.fl, ep_lat, ep_lon, loc_grid(1), loc_grid(2), loc_sta(1), loc_sta(2), rr, tt);
            timeshift(count) = output_temp.traveltime_actual; % -timeshift0;
            %timeshift_temp;
            timeshift_flat(count) = output_temp.traveltime_flat;

            %new_grid_x(i, j) = loc_new_grid(2);
            %new_grid_y(i, j) = loc_new_grid(1)
            
            %new_grid_timeshift(i, j, ii) = timeshift(count);
            %old_grid_timeshift(i, j, ii) = timeshift_flat(count);
            count = count + 1;
        end
    end
    [gridx, gridy] = ndgrid(bpgrid_x_vec, bpgrid_y_vec);
    F1 = scatteredInterpolant(loc_new_grid_x(:), loc_new_grid_y(:), timeshift(:));
    F2 = scatteredInterpolant(loc_flat_grid_x(:), loc_flat_grid_y(:), timeshift_flat(:));

    v1 = F1(gridx, gridy);
    v2 = F2(gridx, gridy);    
    ep_newtimeshift = F1(ep_lon, ep_lat);
    ep_oldtimeshift = F2(ep_lon, ep_lat);
    
    ret1.timeshift(:, :, ii) = v1(:, :); %- ep_newtimeshift; % v1(xslices/2+1, yslices/2+1); % new_grid_timeshift(ii, :, :);
    ret1.timeshift_old(:, :, ii) = v2(:, :); %v2(xslices/2+1, yslices/2+1);
    ret1.lat(:, :, ii) = gridy;
    ret1.lon(:, :, ii) = gridx;
    ret1.stalat(ii) = ret.lat(ii);
    ret1.stalon(ii) = ret.lon(ii);
    ret1.ep_newtraveltime(ii) = ep_newtimeshift;
    ret1.ep_flattraveltime(ii) = ep_oldtimeshift;
end

ret1.lat_step = lat_step;
ret1.lon_step = lon_step;

save('timeshift_newmoho_eu2.mat', 'ret1', '-v7.3');



