clear;
[month,day,hour,min,sec,lat,lon]=textread('newtime.txt','2006-%d-%d:%d:%d:%f %f %f');
MonthDays=[31 28 31 30 31 30 31 31 30 31 30 31]
timeinyear=zeros(size(month));
for ii=1:size(month)
    monthdaysum=0; 
    year=2006;  
    for i=1:month(ii)-1
          monthdaysum=monthdaysum+MonthDays(i);
    end

    if mod(year,4)~=0
    lengthofyear=365;
    end
    if mod(year,4)==0&&month(ii)<=2
    lengthofyear=366;

    else  mod(year,4)==0&&month(ii)>2
     lengthofyear=366;
     monthdaysum=monthdaysum+1;
    end
    timeinyear(ii)=year+((monthdaysum+day(ii))+hour(ii)/24+min(ii)/24/60+sec(ii)/24/60/60)/lengthofyear;

end
%timeinyearmore=timeinyear(1):1/24/60/365:timeinyear(size(month));
%latmore=interp1(timeinyear,lat,timeinyearmore)';
%lonmore=interp1(timeinyear,lon,timeinyearmore)';



[s1,s2,s3,s4,s5,s6]=textread('tempuwt','%s %s %s %s %s %s');
for i=1:5869*6
year0(i)=0;
month0(i)=0;
day0(i)=0;
hour0(i)=0;
min0(i)=0;
mag0(i)=0;
end
for ii=1:5869
    ss1='';
    ss1=cell2mat(s1(ii));
    year0(ii*6-5)=str2double(ss1(1:2))+2000;
    month0(ii*6-5)=str2double(ss1(3:4));
    day0(ii*6-5)=str2double(ss1(5:6));
    hour0(ii*6-5)=str2double(ss1(8:9));
    min0(ii*6-5)=str2double(ss1(10:11));
    mag0(ii*6-5)=str2double(ss1(12:18));
    
    
    ss2='';
    ss2=cell2mat(s2(ii));
    if(size(ss2)~=size(''))
    year0(ii*6-4)=str2double(ss1(1:2))+2000;
    month0(ii*6-4)=str2double(ss1(3:4));
    day0(ii*6-4)=str2double(ss1(5:6));
    hour0(ii*6-4)=str2double(ss2(1:2));
    min0(ii*6-4)=str2double(ss2(3:4));
    mag0(ii*6-4)=str2double(ss2(5:11));
    end
    
    ss3='';
    ss3=cell2mat(s3(ii));
    if(size(ss3)~=size(''))
    year0(ii*6-3)=str2double(ss1(1:2))+2000;
    month0(ii*6-3)=str2double(ss1(3:4));
    day0(ii*6-3)=str2double(ss1(5:6));
    hour0(ii*6-3)=str2double(ss3(1:2));
    min0(ii*6-3)=str2double(ss3(3:4));
    mag0(ii*6-3)=str2double(ss3(5:11));
    end
    ss4='';
    ss4=cell2mat(s4(ii));
    if(size(ss4)~=size(''))
    year0(ii*6-2)=str2double(ss1(1:2))+2000;
    month0(ii*6-2)=str2double(ss1(3:4));
    day0(ii*6-2)=str2double(ss1(5:6));
    hour0(ii*6-2)=str2double(ss4(1:2));
    min0(ii*6-2)=str2double(ss4(3:4));
    mag0(ii*6-2)=str2double(ss4(5:11));
    end
    
    ss5='';
    ss5=cell2mat(s5(ii));
    if(size(ss5)~=size(''))
    year0(ii*6-1)=str2double(ss1(1:2))+2000;
    month0(ii*6-1)=str2double(ss1(3:4));
    day0(ii*6-1)=str2double(ss1(5:6));
    hour0(ii*6-1)=str2double(ss5(1:2));
    min0(ii*6-1)=str2double(ss5(3:4));
    mag0(ii*6-1)=str2double(ss5(5:11));
    end
    
    ss6='';
    ss6=cell2mat(s6(ii));
    if(size(ss6)~=size(''))
    year0(ii*6)=str2double(ss1(1:2))+2000;
    month0(ii*6)=str2double(ss1(3:4));
    day0(ii*6)=str2double(ss1(5:6));
    hour0(ii*6)=str2double(ss6(1:2));
    min0(ii*6)=str2double(ss6(3:4));
    mag0(ii*6)=str2double(ss6(5:11));
    end
end
for ii=1:5869*6
    monthdaysum=0; 
    year=2006;  
    for i=1:month0(ii)-1
          monthdaysum=monthdaysum+MonthDays(i);
    end

    if mod(year,4)~=0
    lengthofyear=365;
    end
    if mod(year,4)==0&&month0(ii)<=2
    lengthofyear=366;

    else  mod(year,4)==0&&month0(ii)>2
     lengthofyear=366;
     monthdaysum=monthdaysum+1;
    end
    timeinyear0(ii)=year+((monthdaysum+day0(ii))+hour0(ii)/24+min0(ii)/24/60)/lengthofyear;

end
%timeinyearmore=timeinyear(1):1/24/60/365:timeinyear(size(month));
lat0=interp1(timeinyear,lat,timeinyear0)';
lon0=interp1(timeinyear,lon,timeinyear0)';
abc=[timeinyear0' lat0 lon0 mag0']
fid = fopen('abcd.txt','wt');
for ii=1:5869*6
        fprintf(fid,'%f %f %f %f\n',timeinyear0(ii), lat0(ii), lon0(ii) ,mag0(ii));
end
        fclose(fid);
    
    
    
    
