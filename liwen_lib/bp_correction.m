function bp_correction(evt_lat, evt_lon, sta_lat, sta_lon)
%this function should have input the outgoing ray direction vector and
%point, then backtrace it and get its direction and point at some depth
%then use this point and depth to compute the actual trace using realistic
%geometry beneath the surface. to get the time shift with respect to the
%original one.

     dep_flat = ;
     lat_flat = ;
     lon_flat = ;
     
     ret = taupPierce('prem', 20, 'P', 'evt', [0, 0], 'sta', [45, 45], 'pierce', 200, 'nodiscon')
     lat_ori = ret.pierce.latitude;
     lon_ori = ret.pierce.longitude;
     
end