%function plot_three_layers(grid_lat_range, grid_lon_range, v0_point, moho_norm_vec, third_layer_depth)
    grid_lat_range = [-1, 1];
    grid_lon_range = [-1, 1];
    v0_point = [0, 0, -50];
    moho_norm_vec = [0.1729,   0.0806,    0.9816];
    d = -dot(moho_norm_vec, v0_point);
    third_layer_depth = 80;
    ini_point = [0,0,-20];
    
    d2km = 6400*2*3.1415926/360;
    km2d = 1 / d2km;

    range_lat = (grid_lat_range(2) - grid_lat_range(1))*d2km;
    range_lon = (grid_lon_range(2) - grid_lon_range(1))*d2km;
    interval_lat = range_lat/100;
    interval_lon = range_lon/100;
    %xx = linspace();
    [xx, yy] = ndgrid(grid_lat_range(1)*d2km:interval_lat:grid_lat_range(2)*d2km, grid_lon_range(1)*d2km:interval_lon:grid_lon_range(2)*d2km);

    third_layer_zz = -third_layer_depth * ones(size(xx, 1), size(xx, 2));
    moho_zz = (-moho_norm_vec(1)*xx - moho_norm_vec(2)*yy - d) / moho_norm_vec(3);
    moho_flat_zz = v0_point(3) * ones(size(xx,1), size(xx, 2));
    ini_zz = ini_point(3) * ones(size(xx, 1), size(yy, 2));
    figure;
    hold on;
    surf(xx*km2d, yy*km2d, third_layer_zz);
    surf(xx*km2d, yy*km2d, moho_zz);
    surf(xx*km2d, yy*km2d, moho_flat_zz);
    surf(xx*km2d, yy*km2d, ini_zz);

%     ini_ray_x = linspace(ini_point(1), I0_0(1), 50);
%     ini_ray_y = linspace(ini_point(2), I0_0(2), 50);
%     ini_ray_tan = ini_vec(3) / sqrt(ini_vec(1)^2 + ini_vec(2)^2);
%     ini_ray_z = ini_ray_tan * sqrt((ini_ray_x-ini_point(1)).^2 + (ini_ray_y-ini_point(2)).^2);
%     plot3(ini_ray_x, ini_ray_y, ini_ray_z+ini_point(3));
% 
% 
%     refr_ray_x = linspace(I0_0(1), I0_1(1), 50);
%     refr_ray_y = linspace(I0_0(2), I0_1(2), 50);
%     refr_ray_tan = refracted_vec0(3) / sqrt(refracted_vec0(1)^2 + refracted_vec0(2)^2);
%     refr_ray_z = refr_ray_tan * sqrt((refr_ray_x-I0_0(1)).^2 + (refr_ray_y-I0_0(2)).^2);
%     plot3(refr_ray_x, refr_ray_y, refr_ray_z+I0_0(3));
%end