% plot the new timeshift
close all;
figure;
hold all;
scatter(loc_new_grid_x, loc_new_grid_y, 20);
scatter(loc_old_grid_x, loc_old_grid_y, 30);

[gridx, gridy] = ndgrid(ep_lon+bparea_span(1,1):lon_step:ep_lon+bparea_span(1,2), ep_lat+bparea_span(2,1):lat_step:ep_lat+bparea_span(2,2));
F1 = scatteredInterpolant(loc_new_grid_x(:), loc_new_grid_y(:), timeshift(:));
F2 = scatteredInterpolant(loc_old_grid_x(:), loc_old_grid_y(:), timeshift_old(:));

v1 = F1(gridx, gridy);
v2 = F2(gridx, gridy);
figure;

surf(gridx, gridy, v1);
colorbar;
title('new\_timeshift');
figure;
surf(gridx, gridy, v2);
title('old\_timeshift');
colorbar;