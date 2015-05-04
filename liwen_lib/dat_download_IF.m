function dat_download_IF(yr, mo, dy, hr, mi, se, duration)
   % dat_download_IF is used to download data from DMC using irisFetch 
    javaaddpath('./IRIS-WS-2.0.12.jar');
    %basetime_num = datenum(str2num(yr), str2num(mo), str2num(dy),str2num(hr),...
    %    str2num(mi), str2num(se));
    
    basetime_num = datenum(yr, mo, dy, hr, mi, se);
    npieces = ceil(duration/120);
    
    for i = 1:npieces
        bt_str = datestr(basetime_num+(i-1)*1/12, 31);
        
        if i == npieces
            et_str = datestr(basetime_num+duration/120*1/12, 31);
        else
            et_str = datestr(basetime_num+i*1/12, 31);
        end
        stations = irisFetch.Stations('station', 'TA', '', '*', 'BHZ', 'starttime',...
            bt_str, 'endtime', et_str, 'matchtimeseries', true);
        nsta = numel(stations);
    
        eff_i = 0;
        clear Traces;
        for idx = 1:nsta
            sta = stations(idx);
            display(sta.StationCode);
            trace = irisFetch.Traces(sta.NetworkCode, sta.StationCode, '*',...
                'BHZ', bt_str, et_str);
        
            if isempty(trace)
                display([sta.StationCode, '', num2str(i), ' is empty!!']);
            else
                eff_i = eff_i + 1;
                Traces(eff_i) = trace(1);  % sometime trace would contain 2 structs so that is why I do this
                size(trace)
                size(Traces(eff_i))
            end
        end
        bt_datestr1 = datestr(basetime_num+(i-1)*1/12, 30);
        if i == npieces
            et_datestr1 = datestr(basetime_num+duration/120*1/12, 30);
        else
            et_datestr1 = datestr(basetime_num+i*1/12, 30);
        end
        file_nm = ['IF', bt_datestr1, '_', et_datestr1, '.mat'];
        save(file_nm, 'Traces', '-v7.3');
    end
    display('save done!');
end