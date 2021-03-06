function [I1_2, refracted_vec, I0_2] = raytracing_liwen_upward(ini_point, ini_vec, v0_point, moho_norm_vec, evt_point)
    %ini_point = ini_point_in;
    %v0_point = [x_v0, y_v0, -dep_v0];
    %moho_norm_vec = [x_norm, y_norm, z_norm];
    %ini_vec = [ini_vec_x, ini_vec_y, ini_vec_z];
    P1_point = ini_point + ini_vec;
    
    n_index = 8.6/6.2;  
    moho_norm_flat_vec = [0, 0, 1];
    [I, check] = plane_line_intersect(moho_norm_vec, v0_point, ini_point, P1_point);
    [I0, check0] = plane_line_intersect(moho_norm_flat_vec, v0_point, ini_point, P1_point);

    refracted_vec = refracted(ini_vec, moho_norm_vec, n_index);
    refracted_vec0 = refracted(ini_vec, [0, 0, 1], n_index);

    surface_point = evt_point;% [0, 0, -evt_depth];
    surface_norm_vec = [0, 0, 1];
    
    I1_2 = plane_line_intersect(surface_norm_vec, surface_point, I,  refracted_vec+I);
    I0_2 = plane_line_intersect(surface_norm_vec, surface_point, I0, refracted_vec0+I0);

    d = -dot(moho_norm_vec, v0_point);
    [xx, yy] = ndgrid(-50:50, -50:50);
    
    surface_zz = surface_point(3) * ones(size(xx,1), size(xx, 2)); %0*xx + 0*yy;
    moho_zz = (-moho_norm_vec(1)*xx - moho_norm_vec(2)*yy - d) / moho_norm_vec(3);
    moho_flat_zz = v0_point(3)*ones(size(xx,1), size(xx, 2));
    inci_zz = ini_point(3) * ones(size(xx, 1), size(yy, 2));
    figure;
    hold on;
    surf(xx, yy, surface_zz);
    surf(xx, yy, moho_zz);
    surf(xx, yy, moho_flat_zz);
    surf(xx, yy, inci_zz);
    
    %[ini_ray_x, ini_ray_y] = ndgrid(ini_point(1):.1:I(1), ini_point(2):.1:I(2));
    ini_ray_x = linspace(ini_point(1), I(1), 50);
    ini_ray_y = linspace(ini_point(2), I(2), 50);
    ini_ray_tan = ini_vec(3) / sqrt(ini_vec(1)^2 + ini_vec(2)^2);
    ini_ray_z = ini_ray_tan * sqrt((ini_ray_x-ini_point(1)).^2 + (ini_ray_y-ini_point(2)).^2);
    plot3(ini_ray_x, ini_ray_y, ini_ray_z+ini_point(3));
       
    ini2_ray_x = linspace(ini_point(1), I0(1), 50);
    ini2_ray_y = linspace(ini_point(2), I0(2), 50);
    ini2_ray_tan = ini_vec(3) / sqrt(ini_vec(1)^2 + ini_vec(2)^2);
    ini2_ray_z = ini_ray_tan * sqrt((ini2_ray_x-ini_point(1)).^2 + (ini2_ray_y-ini_point(2)).^2);
    plot3(ini2_ray_x, ini2_ray_y, ini2_ray_z+ini_point(3));
    
    %
    
    refr_ray_x = linspace(I(1), I1_2(1), 50);
    refr_ray_y = linspace(I(2), I1_2(2), 50);
    refr_ray_tan = refracted_vec(3) / sqrt(refracted_vec(1)^2 + refracted_vec(2)^2);
    refr_ray_z = refr_ray_tan * sqrt((refr_ray_x-I(1)).^2 + (refr_ray_y-I(2)).^2);
    plot3(refr_ray_x, refr_ray_y, refr_ray_z+I(3));
    
    
    refr2_ray_x = linspace(I0(1), I0_2(1), 50);
    refr2_ray_y = linspace(I0(2), I0_2(2), 50);
    refr2_ray_tan = refracted_vec0(3) / sqrt(refracted_vec0(1)^2 + refracted_vec0(2)^2);
    refr2_ray_z = refr2_ray_tan * sqrt((refr2_ray_x-I0(1)).^2 + (refr2_ray_y-I0(2)).^2);
    plot3(refr2_ray_x, refr2_ray_y, refr2_ray_z+I0(3));
    
    
    shift = I0_2 - I1_2
end