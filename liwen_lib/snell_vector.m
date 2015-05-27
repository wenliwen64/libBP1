function refract_vec = snell_vector( inci_vec, normal_vec, relative_n)
    %relative_n = v_ori/v_next;
    refract_vec = relative_n * (cross(normal_vec, cross(-normal_vec, inci_vec))) - normal_vec*(sqrt(1 - (relative_n*relative_n) * dot(cross(normal_vec, inci_vec), cross(normal_vec, inci_vec))));
end