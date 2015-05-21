function [loc_third_point, refracted_vec0, timeshift] = raytracing_liwen(ini_point, ini_vec, v0_point, moho_norm_vec, third_layer_depth, nindex, varargin)
% Calculate and plot the ray-tracing
% Input: 
%        ini_point: take off point of the initial ray;
%        ini_vec: initial ray's directional vector;
%        v0_point: one point on the middle layer;
%        moho_norm_vec: norm vec on the v0_point;
%        third_layer_depth: third layer depth/ the receiving layer;
%        nindex: the n_index on the refracting surface;
% Output: 
%        loc_third_point: the location of the end point of the ray;
%        refracted_vec0: refracted ray's directional vector;
% Example:
%        loc_third_point, refracted_vec0] = raytracing_liwen_downward([0,0,0], [1,1,-6], [0,0,-20], [0,0,1], 80, 6.2/8.6, 'direction', 'up', 'pic', 'yes');
    d2km = 6400*2*3.1415926/360;
    km2d = 1 / d2km;
    P1_point = ini_point + ini_vec;
    
    n_index = nindex; % 6.2 / 8.6;
    moho_norm_vec_flat = [0, 0, 1];
    % Calculate the intersection point
    [I0_0, check0] = plane_line_intersect(moho_norm_vec, v0_point, ini_point, P1_point);
    [I0_0_flat, check_flat] = plane_line_intersect(moho_norm_vec_flat, v0_point, ini_point, P1_point);
    % Calculate the refracted ray's directional vector
    refracted_vec0 = refracted(ini_vec, moho_norm_vec, n_index);
    refracted_vec0_flat = refracted(ini_vec, moho_norm_vec_flat, n_index);

    third_layer_point = [0, 0, -third_layer_depth];
    third_layer_point_flat = [0, 0, -third_layer_depth];
    third_layer_norm_vec = [0, 0, 1];
    third_layer_norm_vec_flat = [0, 0, 1];
    
    %I1_2 = plane_line_intersect(surface_norm_vec, surface_point, I,  refracted_vec+I);
    I0_1 = plane_line_intersect(third_layer_norm_vec, third_layer_point, I0_0, refracted_vec0 + I0_0);
    I0_1_flat = plane_line_intersect(third_layer_norm_vec, third_layer_point_flat, I0_0_flat, ...
        refracted_vec0_flat + I0_0_flat);
    loc_third_point = [I0_1(1), I0_1(2), -third_layer_depth];
    loc_third_point_flat = [I0_1_flat(1), I0_1_flat(2), -third_layer_depth];

    % Calculate the timeshift;
    l1 = norm(ini_point - I0_0);
    l2 = norm(I0_1 - I0_0);
    l1_flat = norm(ini_point - I0_0_flat);
    l2_flat = norm(I0_1_flat - I0_0_flat);
    dl1 = l1 - l1_flat;
    dl2 = l2 - l2_flat;
    
    if(strcmp(varargin{1}, 'direction') && strcmp(varargin{2}, 'down'))
        v1 = 6.2;
        v2 = 8.6;
    elseif(strcmp(varargin{1}, 'direction') && strcmp(varargin{2}, 'up'))
        v1 = 8.6;
        v2 = 6.2;
    else 
        error('bad direction!');
    end
    
    timeshift = dl1 / v1 + dl2 /v2;
    
    if nargin <= 6 
        return;
    end
    
    ii = 3;
    while ii<=length(varargin)
        switch lower(varargin{ii})
            case {'pic'}
                if ~(isa(varargin{ii+1}', 'char'))
                    error('Input string ''yes'' or ''no'' as option of drawing picture!')
                end
                
                if varargin{ii+1} == 'yes'
                    d = -dot(moho_norm_vec, v0_point);
                    [xx, yy] = ndgrid(-50:50, -50:50);

                    third_layer_zz = -third_layer_depth * ones(size(xx, 1), size(xx, 2));
                    moho_zz = (-moho_norm_vec(1)*xx - moho_norm_vec(2)*yy - d) / moho_norm_vec(3);
                    moho_flat_zz = v0_point(3) * ones(size(xx,1), size(xx, 2));
                    ini_zz = ini_point(3) * ones(size(xx, 1), size(yy, 2));
                    figure;
                    hold on;
                    surf(xx, yy, third_layer_zz);
                    surf(xx, yy, moho_zz);
                    surf(xx, yy, moho_flat_zz);
                    surf(xx, yy, ini_zz);

                    ini_ray_x = linspace(ini_point(1), I0_0(1), 50);
                    ini_ray_y = linspace(ini_point(2), I0_0(2), 50);
                    ini_ray_tan = ini_vec(3) / sqrt(ini_vec(1)^2 + ini_vec(2)^2);
                    ini_ray_z = ini_ray_tan * sqrt((ini_ray_x-ini_point(1)).^2 + (ini_ray_y-ini_point(2)).^2);
                    plot3(ini_ray_x, ini_ray_y, ini_ray_z+ini_point(3));


                    refr_ray_x = linspace(I0_0(1), I0_1(1), 50);
                    refr_ray_y = linspace(I0_0(2), I0_1(2), 50);
                    refr_ray_tan = refracted_vec0(3) / sqrt(refracted_vec0(1)^2 + refracted_vec0(2)^2);
                    refr_ray_z = refr_ray_tan * sqrt((refr_ray_x-I0_0(1)).^2 + (refr_ray_y-I0_0(2)).^2);
                    plot3(refr_ray_x, refr_ray_y, refr_ray_z+I0_0(3));
                else
                end
                ii = ii+2;
            otherwise
                error('Unknow Option %s \n', varargin{ii});
        end
    end
    
        
    %shift = I0_2 - I1_2
end