function refracted_vec = refracted(ini_vec, normal_vec, n_index)

    cos_i = dot(ini_vec, normal_vec) / norm(ini_vec) / norm(normal_vec);

    
    if cos_i <= 0
        normal_vec = -normal_vec;    
    end
     
    sin_i = sqrt(1 - cos_i^2);
    
    sin_r = 1/n_index * sin_i;
    cos_r = sqrt(1 - sin_r^2);
    
    cross_vec = cross(ini_vec, normal_vec);
    
    if cross_vec == 0
        refracted_vec = ini_vec / norm(ini_vec);
    else
        surface_vec = cross(normal_vec, cross_vec);
      
        refracted_x = 1*cos_r*normal_vec(1)/norm(normal_vec) + 1*sin_r*surface_vec(1)/norm(surface_vec);
        refracted_y = 1*cos_r*normal_vec(2)/norm(normal_vec) + 1*sin_r*surface_vec(2)/norm(surface_vec);
        refracted_z = 1*cos_r*normal_vec(3)/norm(normal_vec) + 1*sin_r*surface_vec(3)/norm(surface_vec);
    
        refracted_vec = [refracted_x, refracted_y, refracted_z];
    end
end