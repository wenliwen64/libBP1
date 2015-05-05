function usarray_mat2bp0(filename, lat_ep, lon_ep, dep_ep)
% 38.297,  142.372 is the tohoku earthquake location
    display(filename);
    load(filename);
    %calibration event location
    %this is for nepal_2015 main shock
    lon0 = lon_ep;
    lat0 = lat_ep;
    dep0 = dep_ep;
    
    nsta = numel(Traces);
    npts = numel(Traces(1).data);

    r = nan(nsta ,2);
    lat = nan(1, nsta);
    lon = nan(1, nsta);
    rdis = nan(1, nsta);
    az = nan(1, nsta);
    sta_nm = char(zeros(nsta, 8));
    t1 = nan(1, nsta);
    
    npts_dec = numel(decimate(Traces(1).data, 4));
    npts_dec
    xori = zeros(nsta, npts_dec);
    size(xori)
    
    for j = 1:nsta
        lon(j) = Traces(j).longitude;
        r(j, 1) = lon(j);
        lat(j) = Traces(j).latitude;
        r(j, 2) = lat(j);
        
        len_nm = numel(Traces(j).station);
        sta_nm(j, 1:len_nm) = Traces(j).station;
        
        down_sample = decimate(Traces(j).data, 4);
        down_sample = down_sample - mean(down_sample);
        npts_temp = numel(down_sample);
        xori(j, 1:npts_temp) = down_sample;
  
        t1(j) = Traces(j).startTime;
        [rdis(j), az(j)] = distance(r(j, 2), r(j, 1), lat0, lon0);
    end
    
    new_str.r = r;
    new_str.nm = sta_nm;
    new_str.lat = lat;
    new_str.lon = lon;
    new_str.rdis = rdis;
    new_str.az = az;
    new_str.t1 = t1;
    new_str.sr = 10;
    new_str.n = nsta;
    new_str.xori = xori;
    new_str.lat0 = lat0;
    new_str.lon0 = lon0;
    new_str.dep0 = dep0;
    
    [f_dir, f_nm, f_ext] = fileparts(filename);
    save([f_nm '.v1.mat'], 'new_str');
end