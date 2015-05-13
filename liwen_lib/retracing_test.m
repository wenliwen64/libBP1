function loc_new_grid = retracing_test(loc_sta, loc_grid, depth)
     ret = taupTime('prem', depth, 'P', 'sta', loc_sta, 'evt', loc_grid);
     sin_takeoff_angle = ret.rayParam/(6400*2*3.1415926/360)*6.2/57.3;
     cos_takeoff_angle = sqrt(1 - sin_takeoff_angle^2);
     [dk, dd, daze, dazs] = distaz(loc_sta(1), loc_sta(2), loc_grid(1), loc_grid(2));
     ini_vec = [1*sin_takeoff_angle*sin(daze), 1*sin_takeoff_angle*cos(daze), -1*cos_takeoff_angle];
     %ret2 = taupPierce('prem', depth, 'P', 'sta', loc_sta, 'evt', loc_grid);
     ini_point = [loc_grid(1), loc_grid(2), -20];
     v0_point = [loc_grid(1), loc_grid(2), -40];
     moho_norm_flat_vec = [0, 0, -1];
     third_layer_depth = 80;
     moho_norm_actual_vec = [1, 1, 6];
     [loc_third_point, downward_vec] = raytracing_liwen_downward(ini_point, ini_vec, v0_point, ...
         moho_norm_flat_vec, third_layer_depth);
     [loc_new_grid, new_vec, I0_2] = raytracing_liwen_upward(loc_third_point, -downward_vec, v0_point, moho_norm_actual_vec, [loc_grid(1), loc_grid(2), -depth]);
     ret.rayParam
     sin_takeoff_angle
     cos_takeoff_angle
     I0_2
end