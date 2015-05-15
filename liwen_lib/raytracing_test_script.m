%raytracing functions test script
ini_point = [15, 75,-15];
ini_vec = [1, 0, -5];
v0_point = [20,80, -50];
moho_norm_vec_flat = [0, 0, 1];
moho_norm_actual_vec = [1, 0, 6];

nindex_down = 6.2/8.6;
nindex_up = 1 / nindex_down;
third_layer_depth = 80;
[loc_third_point, refracted_vec0, timeshift] = raytracing_liwen_downward(ini_point,...
    ini_vec, v0_point, moho_norm_vec_flat, third_layer_depth, nindex_down, 'direction', 'down');
[loc_new_grid, upward_vec, timeshift_up] = raytracing_liwen_downward(loc_third_point, -refracted_vec0, ...
         v0_point, moho_norm_actual_vec, 15, nindex_up, 'direction', 'up');