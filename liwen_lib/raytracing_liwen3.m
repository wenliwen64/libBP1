function output = raytracing_liwen3(input)
% Calculate and plot the ray-tracing
% Input(r-[lon, lat] in km, loc-[lat, lon] in degree): 
%        input.ini_point_r: take off point of the initial ray;
%        input.ini_vec: initial ray's directional vector, in km;
%        input.v0_point_r: one point on the middle layer;
%        input.moho_norm_vec: norm vec on the v0_point, the actual one, the one we want to test, in km;
%        input.third_layer_depth_retracing: third layer depth/ the receiving layer;
%        input.epicenter_r: epicenter location;
%        input.bprange_r: back-projection area;
%        //input.nindex: the n_index on the refracting surface;
%        input.direction: 1-upward, 0-downward, upward only here;
%        input.plot: 1-plot, 0-not.
% Output: 
%        output.third_point_r: the [x, y](in km, lon = x, lat =y) coordinates;
%        output.third_point_flat_r: the [x, y] coordinates for third point
%            using flat moho surface(in km);
%        output.refracted_vec: refracted ray's directional vector using actual moho surface(in km);
%        output.refracted_vec_flat: refracted ray's directional vector using flat moho surface(in km);
%        output.traveltime: travel time for actual moho surface;
%        output.traveltime_flat: travel time for flat moho surface;
% Example:
%        input.ini_point_r = [0,0,-80];
%        input.ini_vec = [1,1,6];
%        input.v0_point_r = [0,0,-50];
%        input.moho_norm_vec = [0.1729, 0.0806, 0.9816];
%        input.third_layer_depth_retracing = 20;
%        input.epicenter_r = [0,0];
%        input.bprange_r = [-1,1;-1,1];
%        input.direction = 1;
%        input.plot = 0;
%        output = raytracing_liwen(input);
%        //output = raytracing_liwen([0,0,0], [1,1,-6], [0,0,-20], [0,0,1], 80, 6.2/8.6, [0,0], [-1,1;-1,1],'direction', 'up', 'pic', 'yes');
    ini_point_r = input.ini_point_r;
    ini_vec = input.ini_vec;
    v0_point_r = input.v0_point_r;
    moho_norm_vec = input.moho_norm_vec;
    third_layer_depth = input.third_layer_depth_retracing;
    epicenter_r = input.epicenter_r;
    bprange_r = input.bprange_r;
    
    d2km = 6371*2*3.1415926/360;
    km2d = 1 / d2km;
    P1_point = ini_point_r + ini_vec;
    
    %n_index = nindex; % 6.2 / 8.6;
    v_up = 6.8;
    v_bot = 8.2;
    moho_norm_vec_flat = [0, 0, -1];
    
    % Calculate the intersection point
    [I0_0, check0] = plane_line_intersect(moho_norm_vec, v0_point_r, ini_point_r, P1_point);
    [I0_0_flat, check_flat] = plane_line_intersect(moho_norm_vec_flat, v0_point_r, ini_point_r, P1_point);
    
    % Calculate the refracted ray's directional vector using Liwen's
    % refraction code, the third parameter should be the n_index(> 0), so
    % here we only ALLOW upward.
    refracted_vec0 = refracted(ini_vec, moho_norm_vec, v_bot/v_up);
    refracted_vec0_flat = refracted(ini_vec, moho_norm_vec_flat, v_bot/v_up);

    third_layer_point = [0, 0, -third_layer_depth];
    third_layer_norm_vec = [0, 0, 1]; % the earthquake layer
    
    %I1_2 = plane_line_intersect(surface_norm_vec, surface_point, I,  refracted_vec+I);
    I0_1 = plane_line_intersect(third_layer_norm_vec, third_layer_point, I0_0, I0_0+refracted_vec0);
    I0_1_flat = plane_line_intersect(third_layer_norm_vec, third_layer_point, I0_0_flat, I0_0_flat+refracted_vec0_flat );
    
    third_point_r = [I0_1(1), I0_1(2), -third_layer_depth];
    third_point_flat_r = [I0_1_flat(1), I0_1_flat(2), -third_layer_depth];

    % Calculate the traveltime;
    l1 = norm(ini_point_r - I0_0);
    l2 = norm(I0_1 - I0_0);
    l1_flat = norm(ini_point_r - I0_0_flat);
    l2_flat = norm(I0_1_flat - I0_0_flat);
    
    if(input.direction == 0) % 0 = downward
        v1 = 6.8;
        v2 = 8.2;
    elseif(input.direction == 1) % 1 = upward
        v1 = 8.2;
        v2 = 6.8;
    else 
        error('bad direction!');
    end
    
    traveltime = l1 / v1 + l2 /v2;
    traveltime_flat = l1_flat / v1 + l2_flat / v2;
    
    output.third_point_r = third_point_r;
    output.third_point_flat_r = third_point_flat_r;
    output.refracted_vec = refracted_vec0;
    output.refracted_vec_flat = refracted_vec0_flat;
    output.traveltime = traveltime;
    output.traveltime_flat = traveltime_flat;
                 
    if input.plot == 1
        lon_step = (bprange_r(1,2) - bprange_r(1,1)) / 40;
        lat_step = (bprange_r(2,2) - bprange_r(2,1)) / 40;
        [xx, yy] = ndgrid(epicenter_r(1)+bprange_r(1,1):lon_step:epicenter_r(1)+bprange_r(1,2), ...
            epicenter_r(2)+bprange_r(2,1):lat_step:epicenter_r(2)+bprange_r(2,2)); % ndgrid(-50:50, -50:50);
        
        third_layer_zz = -third_layer_depth * ones(size(xx, 1), size(xx, 2));
        
        d = -dot(moho_norm_vec, v0_point_r);
        moho_zz = (-moho_norm_vec(1)*xx - moho_norm_vec(2)*yy - d) / moho_norm_vec(3);
        moho_flat_zz = v0_point_r(3) * ones(size(xx,1), size(xx, 2));
        ini_zz = ini_point_r(3) * ones(size(xx, 1), size(yy, 2));
        
        %figure;
        hold on;
        surf(xx, yy, third_layer_zz);
        surf(xx, yy, moho_zz);
        %surf(xx, yy, moho_flat_zz);
        surf(xx, yy, ini_zz);

        ini_ray_x = linspace(ini_point_r(1), I0_0(1), 50);
        ini_ray_y = linspace(ini_point_r(2), I0_0(2), 50);
        ini_ray_tan = ini_vec(3) / sqrt(ini_vec(1)^2 + ini_vec(2)^2);
        ini_ray_z = ini_ray_tan * sqrt((ini_ray_x-ini_point_r(1)).^2 + (ini_ray_y-ini_point_r(2)).^2);
        plot3(ini_ray_x, ini_ray_y, ini_ray_z+ini_point_r(3), 'r*');


        refr_ray_x = linspace(I0_0(1), I0_1(1), 50);
        refr_ray_y = linspace(I0_0(2), I0_1(2), 50);
        refr_ray_tan = refracted_vec0(3) / sqrt(refracted_vec0(1)^2 + refracted_vec0(2)^2);
        refr_ray_z = refr_ray_tan * sqrt((refr_ray_x-I0_0(1)).^2 + (refr_ray_y-I0_0(2)).^2);
        plot3(refr_ray_x, refr_ray_y, refr_ray_z+I0_0(3), 'k');
    else
    end
end