

% function [data,begin,lat,lon,height,nm,sr]=readsm(filename)
function [data begin lat lon nm sr]=readv2(filename,chan,avd)
%chan='UP'
%avd='ACCEL'
file = textread(filename,'%s');
data=NaN;
flag=0;
for i=1:length(file)
    if strcmpi(file(i),'CHAN')&&strcmpi(file(i+2),chan)
        flag=1;
    end
    if (strcmpi(file(i),'TRIGGER')||strcmpi(file(i),'TRIGER')||strcmpi(file(i),'START'))&&flag
        disp([char(file(i+2)) char(file(i+3))])
       begin=datenum([char(file(i+2)) char(file(i+3))]);  
    end
    if strcmpi(file(i),'STATION')&&flag
       nm=file(i+2);
       tmp=char(file(i+3));
       if strcmpi(tmp(end-1),'N')
           sign=1;
       else
           sign=-1;
       end
       lat=sign*str2num(tmp(1:end-2));
       tmp=char(file(i+4));
       if strcmpi(tmp(end),'E')
           sign=1;
       else
           sign=-1;
       end
       lon=sign*str2num(tmp(1:end-1));
    end
    if strcmpi(file(i),'INTERVALS')&&flag
       sr=1/str2num(char(file(i+2)));
    end
    if strcmpi(file(i),'POINTS')&&strcmpi(file(i+2),avd)&&flag
        disp([chan ' ' avd]);
        data=str2num(char(file(i+11)));
        for j=i+12:i+11+str2num(char(file(i-1)))-1
            data=[data str2num(char(file(j)))];
        end
        flag=0;
    end
end
