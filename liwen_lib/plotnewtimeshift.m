% plot the new timeshift
close all;
figure;
hold all;
bpgrid_x_vec = linspace(ep_lon + bparea_span(1, 1), ep_lon + bparea_span(1,2), xslices);
bpgrid_y_vec = linspace(ep_lat + bparea_span(2, 1), ep_lat + bparea_span(2, 2), yslices);
count_temp = 1;
for i = 1:40*40
    
    ret_temp = taupTime('prem', 20, 'P', 'sta', loc_sta, 'evt', [loc_old_grid_y(i), loc_old_grid_x(i)]);
    a(i) = ret_temp.time;
    
end
scatter(loc_old_grid_x, loc_old_grid_y, 10);
scatter(loc_new_grid_x, loc_new_grid_y, 20);
scatter(loc_flat_grid_x, loc_flat_grid_y, 30);

[gridx, gridy] = ndgrid(ep_lon+bparea_span(1,1):lon_step:ep_lon+bparea_span(1,2), ep_lat+bparea_span(2,1):lat_step:ep_lat+bparea_span(2,2));
F1 = scatteredInterpolant(loc_new_grid_x(:), loc_new_grid_y(:), timeshift(:));
F2 = scatteredInterpolant(loc_flat_grid_x(:), loc_flat_grid_y(:), timeshift_old(:));
F3 = scatteredInterpolant(loc_old_grid_x(:), loc_old_grid_y(:), a(:));

v1 = F1(gridx, gridy);
v2 = F2(gridx, gridy);
v3 = F3(gridx, gridy);
figure;

%surf(gridx, gridy, v1-ep_newtimeshift);
surf(gridx, gridy, v1);
colorbar;
title('new\_timeshift');
figure;
%surf(gridx, gridy, v2-ep_oldtimeshift);
surf(gridx, gridy, v2);
title('old\_timeshift');
colorbar;
figure;
surf(gridx, gridy, v1-v2-ep_newtimeshift+ep_oldtimeshift);
%surf(gridx, gridy, v1-v2);
colorbar;
figure;
surf(gridx, gridy, v3);