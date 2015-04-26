

function [data,begin,lat,lon,height,nm,sr]=readsm(filename)
file = textread(filename,'%s','delimiter','\n','whitespace','');
tmp=char(file(6));
nm=tmp(14:end);


tmp=char(file(7));
lat=str2num(tmp(14:end));

tmp=char(file(8));
lon=str2num(tmp(14:end));


tmp=char(file(9));
height=str2num(tmp(18:end));


tmp=char(file(10));
begin=datenum(tmp(14:end));


tmp=char(file(11));
tmp1=strtrim((tmp(18:end)));
sr=str2num(tmp1(1:end-2));



data = str2num(char(file(18)));
for i=19:length(file)
    data=[data str2num(char(file(i)))];
end
data=data*2000/8388608;