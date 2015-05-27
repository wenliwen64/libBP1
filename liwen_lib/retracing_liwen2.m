function [loc_new_grid_flat, loc_new_grid, traveltime_flat, traveltime_actual] = retracing_liwen2(loc_sta, loc_grid, evt_depth, epicenter, bprange) % (lat, lon) tuple
% This function is used to recompute the travel time using different moho
% geomety setup.
% Input: 
%     loc_sta = [lat_sta, lon_sta]
%     loc_grid = [lat_grid, lon_grid]
%     evt_depth = event depth
% Output: 
%     loc_new_grid = location of new grid point after re-raytracing down-up
%                    calculation
%     timeshift = time difference computed by  T_new - T_old
% Example:
%     [loc_new_grid, timeshift] = retracing_liwen([0, 60], [0, 100], 20)

     % nepal mainshock epicenter location
     lat0 = 28.1473; 
     lon0 = 84.7079;
     d2km = 6371*2*3.1415926/360;
     km2d = 1/d2km;
     third_layer_depth = 80;
     ret = taupTime('prem', third_layer_depth, 'P', 'sta', loc_sta, 'evt', loc_grid);
    
     % calculate the initial take off angle
     sin_takeoff_angle = ret(1).rayParam / (6371*2*3.1415926/360) * 8.2 / 57.3; % mysterious constant 57.3, 8.2 is the velocity at 80 depth 
     cos_takeoff_angle = sqrt(1 - sin_takeoff_angle^2);
     [dk, dd, daze, dazs] = distaz(loc_sta(1), loc_sta(2), loc_grid(1), loc_grid(2));
     
     % set up the initial parameters;!!========IMPORTANT=========!!
     ini_vec = [1*sin_takeoff_angle*sind(dazs), 1*sin_takeoff_angle*cosd(dazs), -1*cos_takeoff_angle];% upward
     ini_vec = ini_vec / norm(ini_vec);
     ini_point = [loc_grid(2)*d2km, loc_grid(1)*d2km, -third_layer_depth];% use third_layer as the initial surface
     v0_point = [(lon0)*d2km, (lat0)*d2km, -50];  % point on the moho surface;
     moho_norm_flat_vec = [0, 0, 1];
     third_layer_depth = 80;
     moho_norm_actual_vec = [0.1729, 0.0806, 0.9816]; % based on info at http://earthquake.usgs.gov/earthquakes/eventpage/us20002926#scientific_tensor:us_us_20002926_mwc
     nindex_down = 6.8/8.2; % original_velocity / next_medium_velocity
     nindex_up = 8.2/6.8;
     
     [loc_new_grid_xy_flat, downward_vec, traveltime_flat] = raytracing_liwen2(ini_point, -ini_vec, v0_point, ...
         -moho_norm_flat_vec, evt_depth, nindex_down, epicenter, bprange, 'direction', 'up');
     
     [loc_new_grid_xy, upward_vec, traveltime_actual] = raytracing_liwen(ini_point, -ini_vec, ...
         v0_point, -moho_norm_actual_vec, evt_depth, nindex_up, epicenter, bprange,  'direction', 'up');
     
     loc_new_grid_flat = [loc_new_grid_xy_flat(2)*km2d, loc_new_grid_xy_flat(1)*km2d, evt_depth];% lat, lon
     loc_new_grid = [loc_new_grid_xy(2)*km2d, loc_new_grid_xy(1)*km2d, evt_depth];
     
     traveltime_flat = traveltime_flat + ret.time;
     traveltime_actual = traveltime_actual + ret.time;
     

end