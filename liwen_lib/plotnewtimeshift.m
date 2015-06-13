% plot the new timeshift
%close all;
close all;
figure(100);
hold all;
ep_lat = 28.1473; 
ep_lon = 84.7079;
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
F2 = scatteredInterpolant(loc_flat_grid_x(:), loc_flat_grid_y(:), timeshift_flat(:));
F3 = scatteredInterpolant(loc_old_grid_x(:), loc_old_grid_y(:), a(:));

v1 = F1(gridx, gridy);
v2 = F2(gridx, gridy);
v3 = F3(gridx, gridy);
taup_eptime = F3(ep_lon, ep_lat);
figure(1);

%surf(gridx, gridy, v1-ep_newtimeshift);
surf(gridx, gridy, v1);
colorbar;
title('new\_timeshift');
figure(2);
%surf(gridx, gridy, v2-ep_oldtimeshift);
surf(gridx, gridy, v2);
title('old\_timeshift');
colorbar;

colorbar;
figure(4);
surf(gridx, gridy, v1-v2-ep_newtimeshift+ep_oldtimeshift);
title('variance v1 - v2');
%surf(gridx, gridy, v1-v2);
colorbar;
figure(5);
surf(gridx, gridy, v3);
colorbar;
title('v3');
figure(6);
surf(gridx, gridy, v3-v2-taup_eptime+ep_oldtimeshift);
title('v3 - v2 variance');
colorbar;
figure(7);
surf(gridx, gridy, v3-taup_eptime - v1 +ep_newtimeshift);
title('v3 - v1 variance');
colorbar;
