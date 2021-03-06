function bp_correction(epi_lat, epi_lon, sta_lat, sta_lon, bprange_lat, bprange_lon)
%this function should have input the outgoing ray direction vector and
%point, then backtrace it and get its direction and point at some depth
%then use this point and depth to compute the actual trace using realistic
%geometry beneath the surface. to get the time shift with respect to the
%original one.

     dep_flat = ;
     lat_flat = ;
     lon_flat = ;
     
     nsta = numel(sta_lat);
     
     for i = 1:nsta
           sta_lon_temp = sta_lon(i);
           sta_lat_temp = sta_lat(i);
           loc_sta = [sta_lat_temp sta_lon_temp];
           for j = 1:lat_range
               for k = 1:lon_range
                   loc_grid = [epi_lat + j*dlat, epi_lon + k*dlon];    
                   m_timeshift(i, j, k)  = correct_timeshift(loc_sta, loc_grid, depth, v0_point, moho_norm_vec);
               end
           end
     end
     
     save('m_timeshift.mat', 'm_timeshift');
     ret = taupPierce('prem', 20, 'P', 'evt', [0, 0], 'sta', [45, 45], 'pierce', 200, 'nodiscon')
     lat_ori = ret.pierce.latitude;
     lon_ori = ret.pierce.longitude;     
end