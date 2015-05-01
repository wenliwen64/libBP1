function hinet_rdsac2bp1(filename)
    load(filename);
    opr = readAllSac();
    opr.bpbool = false;
    opr.lon0 = new_str.lon0;
    opr.lat0 = new_str.lat0;
    opr.dep0 = new_str.dep0;
    opr.bp = [0.01 1];
    opr.snrFilterbool = false;
    opr.ori = 300;
    opr.snrFilter = [0.1, 0.5, 2, -20, -10, 100, 130];
    opr.sr = 10;
    
    %new_str has been loaded 
    ret = new_str;
    ret.opr = opr;
    ret.timeshiftall = zeros(1, length(ret.lon));
    ret.lon = round(ret.lon*1e4) / 1e4;
    
    if isfield(ret, 'time')
        ret = rmfield(ret, 'time');
    end
    
    load('ptimes.mat');
    shiftP = interp1(rr, tt, ret.rdis);
    nsta = ret.n;
    
    for i = 1:nsta
        ret.xori(i, :) = specshift(ret.xori(i, :), ret.sr*(shiftP(i) - shiftP(1)));
        ret.nf(i) = std(ret.xori(i, 100:4000));
        ret.xori(i, :) = ret.xori(i, :) / ret.nf(i);
    end
    
    ret.lon0 = ret.opr.lon0;
    ret.lat0 = ret.opr.lat0;
    ret.dep0 = ret.opr.dep0;
    ret.parr = 20;
    ret.begin = 0;
    ret.end = round(size(ret.xori, 2)/ret.sr) - 20;
    ret.step = 1;
    ret.ps = 40;
    ret.qs = 40;
    ret.lonrange = [-3 3];
    ret.latrange = [-4 2];
    ret.fl = 0.5;
    ret.fh = 2;
    ret.fs = 2;
    ret.win = 5;
    ret.dirname = 'nepal_hn';
    ret.Nw = 3;
    
    [f_dir, f_nm, f_ext] = fileparts(filename);
    I = findstr('.', f_nm);
    f_nm = f_nm(1:I(end));
    save([f_nm 'v2.mat'], 'ret', '-v7.3');
end