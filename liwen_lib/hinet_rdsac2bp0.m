function hinet_rdsac2bp0(filename)
    display(filename);
    load(filename);
    %calibration event location
    %this is for nepal_2015 main shock
    lon0 = 84.7079;
    lat0 = 28.1473;
    dep0 = 15;
    
    nsta = numel(sh);
    npts = numel(sh(1).d);

    r = nan(nsta ,2);
    lat = nan(1, nsta);
    lon = nan(1, nsta);
    rdis = nan(1, nsta);
    az = nan(1, nsta);
    sta_nm = char(zeros(nsta, 8));
    t1 = nan(1, nsta);
    
    npts_dec = numel(decimate(sh(1).d, 4));
    npts_dec
    xori = nan(nsta, npts_dec);
    size(xori)
    
    for j = 1:nsta
        lon(j) = sh(j).HEADER.STLO;
        r(j, 1) = lon(j);
        lat(j) = sh(j).HEADER.STLA;
        r(j, 2) = lat(j);
        
        len_nm = numel(sh(j).HEADER.KSTNM);
        sta_nm(j, 1:len_nm) = sh(j).HEADER.KSTNM;
        
        down_sample = decimate(sh(j).d, 4);
        down_sample = down_sample - mean(down_sample);
        npts_temp = numel(down_sample);
        xori(j, 1:npts_temp) = down_sample;
        
        styr = sh(j).HEADER.NZYEAR;
        stdy = sh(j).HEADER.NZJDAY;
        sthr = sh(j).HEADER.NZHOUR;
        stmi = sh(j).HEADER.NZMIN;
        stse = sh(j).HEADER.NZSEC;
   
        dateno = datenum(styr, 1, stdy, sthr, stmi, stse);
        t1(j) = dateno;
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
    save([f_nm '.v0.mat'], 'new_str');    
end
