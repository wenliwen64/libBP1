%test script:

% test the routines.

       figure(100);
            d2km = 6371*2*3.1415926/360;
     km2d = 1/d2km;
       input.ini_point_r = [0,0,-80];
       input.ini_vec = [1,1,6];
       input.v0_point_r = [0,0,-50];
       input.moho_norm_vec = norm_vect_comp(18);
       %input.moho_norm_vec = [0.1729, 0.0806, 0.9816];
       input.third_layer_depth_retracing = 20;
       input.epicenter_r = [0,0];
       input.bprange_r = [-1,1;-1,1]*d2km;
       input.direction = 1;
       input.plot = 1;
       output = raytracing_liwen3(input);