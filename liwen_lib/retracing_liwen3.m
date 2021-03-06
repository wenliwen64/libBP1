function output = retracing_liwen3(input) % (lat, lon) tuple
% This function is used to recompute the travel time using different moho
% geomety setup.
% Input: 
%     input.loc_sta: station's location, [sta_lat, sta_lon] in degree; 
%     input.loc_grid: bp point's location, [lat_grid, lon_grid] in degree;
%     input.evt_depth: event depth in km;
%     input.start_depth: starting point's depth in kilometer;
%     input.third_layer_depth_retracing: the final layer in the wave
%         propagation?
%     input.epicenter_r: epicenter, 2d, [x, y] in km;
%     input.bprange_r: bp area in km, 2d;
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
     
     loc_sta = input.loc_sta;
     loc_grid = input.loc_grid;
     evt_depth = input.evt_depth;
     start_depth = input.start_depth;
     third_layer_depth_retracing = input.third_layer_depth_retracing;
     epicenter_r = input.epicenter_r;
     bprange_r = input.bprange_r;
     moho_norm_flat_vec = input.moho_norm_flat_vec;
     moho_norm_actual_vec = input.moho_norm_actual_vec;%[0.1729, 0.0806, 0.9816];
     
     ret = taupTime('prem', start_depth, 'P', 'sta', loc_sta, 'evt', loc_grid);
    
     % calculate the initial take off angle
     sin_takeoff_angle = ret(1).rayParam / (6371*2*3.1415926/360) * 8.2 / 57.3; % mysterious constant 57.3, 8.2 is the velocity at 80 depth 
     asind(sin_takeoff_angle);
     cos_takeoff_angle = sqrt(1 - sin_takeoff_angle^2);
     [dk, dd, daze, dazs] = distaz(loc_sta(1), loc_sta(2), loc_grid(1), loc_grid(2));
     
     % set up the initial parameters;!!========IMPORTANT=========!!
     ini_vec = [1*sin_takeoff_angle*sind(dazs), 1*sin_takeoff_angle*cosd(dazs), -1*cos_takeoff_angle];% downward
     ini_vec = ini_vec / norm(ini_vec);
     ini_point_r = [loc_grid(2)*d2km, loc_grid(1)*d2km, -start_depth];  % use third_layer as the initial surface
     v0_point_r = [(lon0)*d2km, (lat0)*d2km, -50];  % point on the moho surface;
     %moho_norm_flat_vec = [0, 0, 1];
     %moho_norm_actual_vec = [0.1729, 0.0806, 0.9816];  % based on info at http://earthquake.usgs.gov/earthquakes/eventpage/us20002926#scientific_tensor:us_us_20002926_mw
     
     input_actual.ini_point_r = ini_point_r;
     input_actual.ini_vec = -ini_vec; % upward directing vector
     input_actual.v0_point_r = v0_point_r;
     input_actual.moho_norm_vec = -moho_norm_actual_vec;
     input_actual.third_layer_depth_retracing = third_layer_depth_retracing;
     input_actual.epicenter_r = epicenter_r;
     input_actual.bprange_r = bprange_r;
     input_actual.plot = 0;
     input_actual.direction = 1;
     
     output_actual = raytracing_liwen3(input_actual);
     
     input_flat.ini_point_r = ini_point_r;
     input_flat.ini_vec = -ini_vec;  % upward directing vector
     input_flat.v0_point_r = v0_point_r;
     input_flat.moho_norm_vec = -moho_norm_flat_vec;
     input_flat.third_layer_depth_retracing = third_layer_depth_retracing;
     input_flat.epicenter_r = epicenter_r;
     input_flat.bprange_r = bprange_r;
     input_flat.plot = 0;
     input_flat.direction = 1;
     
     output_flat = raytracing_liwen3(input_flat);
     
     output.loc_new_grid_flat_deg = [output_flat.third_point_r(2)*km2d, output_flat.third_point_r(1)*km2d, evt_depth];% lat, lon
     output.loc_new_grid_deg = [output_actual.third_point_r(2)*km2d, output_actual.third_point_r(1)*km2d, evt_depth];
     
     output.traveltime_flat = output_flat.traveltime  + ret.time;
     output.traveltime_actual = output_actual.traveltime + ret.time;
end