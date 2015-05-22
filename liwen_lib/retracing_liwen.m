function [loc_new_grid, timeshift] = retracing_liwen(loc_sta, loc_grid, evt_depth) % (lat, lon) tuple
     lat0 = 28.1473;%nepal mainshock
     lon0 = 84.7079;
     d2km = 6400*2*3.1415926/360;
     km2d = 1/d2km;
     ret = taupTime('prem', evt_depth, 'P', 'sta', loc_sta, 'evt', loc_grid);
%      ret.rayParam
    
     % calculate the initial take off angle
     sin_takeoff_angle = ret(1).rayParam / (6400*2*3.1415926/360) * 6.8 / 57.3; % mysterious constant 57.3
%      asind(sin_takeoff_angle)
     cos_takeoff_angle = sqrt(1 - sin_takeoff_angle^2);
     [dk, dd, daze, dazs] = distaz(loc_sta(1), loc_sta(2), loc_grid(1), loc_grid(2));
     
     % set up the constant parameters;
     ini_vec = [1*sin_takeoff_angle*sind(dazs), 1*sin_takeoff_angle*cosd(dazs), -1*cos_takeoff_angle];
     ini_point = [loc_grid(2)*d2km, loc_grid(1)*d2km, -evt_depth];
     v0_point = [(lon0)*d2km, (lat0)*d2km, -50];  % moho depth
     moho_norm_flat_vec = [0, 0, -1];
     third_layer_depth = 80;
     
     %moho_norm_actual_vec = [1, 0, 5.1446]; % 11 degree 
     moho_norm_actual_vec = [0.1729,   0.0806,    0.9816];
     nindex_down = 6.2/8.6;
     nindex_up = 8.6/6.2;
     
     [loc_third_point, downward_vec, timeshift_down] = raytracing_liwen(ini_point, ini_vec, v0_point, ...
         moho_norm_flat_vec, third_layer_depth, nindex_down, 'direction', 'down');
     %loc_third_point
     [loc_new_grid, upward_vec, timeshift_up] = raytracing_liwen(loc_third_point, -downward_vec, ...
         v0_point, moho_norm_actual_vec, evt_depth, nindex_up, 'direction', 'up');
     loc_new_grid = [loc_new_grid(2)*km2d, loc_new_grid(1)*km2d, loc_new_grid(3)];
     timeshift = timeshift_down + timeshift_up;
     % checking parameters
%      ret.rayParam
%      downward_vec
%      loc_third_point
%      sin_takeoff_angle
%      cos_takeoff_angle
%      ini_vec
%      daze
%      dazs
%      sin_takeoff_angle
     
%      1*sin_takeoff_angle*cosd(dazs)
end