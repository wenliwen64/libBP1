
function ret=readAllv2(dir,myOpr)
% function to read all sac files in a dir, and output the data and associate information
% usage: function ret=readAllSac(dir,myOpr)
% input: dir the directory where you keep the sac file
%        myOpr is the option object( if there is no input argument the fucntion returns the option object with default settings)
% output: ret is the output object      
% to see the details of the options and the output. read the comment in the
% defination of the these two objects

def = struct(...
        'FileList','filelist',...% the name of the filelist which contains the selected sac files in the dir
         'lat0',0,...% lat of the epicenter
         'lon0',0,...% lon of the epicenter
         'bpbool',false,...% to indicate if you want to bandpass the data for the output 
         'bp',[0.5 2],...% band of the butterworth filter
         'snrFilterbool',false,...% whether to use filtered snr as a select cirteria or not
         'ori',140,...% origin time in the data
         'snrFilter',[0.05 0.5 2 -20 -10 5 10],...% the parameters for snr filter. [1 snr 2 3 band 4 5 time before origin time to define pre signal (noise) 6 7 time after origin for signal]
         'distanceFilterbool',false,...%whether to use distance filter
         'distanceFilter',[26 28.5],...% uper and lower bound distance filter in degrees
         'orderDistancebool',false,...% order the data in distance
         'sr',10,...%sampling rate of the data
         'alignbool',false,...% align the data with initial arrival
         'nrt',0,...% n square root
         'resample',0,...% resample
         'dec',0,...%decimate
         'chan','UP',...%channel
         'avd','ACCEL',...%acceleration, velocity, displacement
         'align',seismoAlign());% options for alignmentCS


    if ~nargin %user wants default options, give it and stop
    ret = def;
    return
    elseif nargin<2% no opr speciified, use default
      Opr=def;
    elseif nargin<3% use the opr specified by user
      Opr=myOpr;
    else 
        error('too many parameters');
    end
%main settings
    FileList=Opr.FileList;
    
    lat0=Opr.lat0;
    lon0=Opr.lon0;
    chan=Opr.chan;
    avd=Opr.avd;
    fid=fopen([dir FileList]);
    count=1;
    r=zeros(1,2);
    nm=zeros(1,4);
    rdis=zeros(1,1);
     az=zeros(1,1);
    t=zeros(1,1);
    x=zeros(1,1);
    t1=zeros(1,1);
    if Opr.resample~=0
   Opr.sr=Opr.sr*Opr.resample;
    end
k=1;
    
% read the sacs with fget_sac    
    while ~feof(fid)
        k=k+1
        
        [B,A]=butter(4,Opr.snrFilter(2:3)/(Opr.sr/2));
        [BB,AA]=butter(4,Opr.bp/(Opr.sr/2));

        dataname=strtrim(fgetl(fid))
        
%         [Zdata begin lats lons height stnm sr ]=readsm([dir dataname]);
         [Zdata begin lats lons stnm sr ]=readv2([dir dataname],chan,avd);
         Zdata=Zdata-median(Zdata);
        
%         [Ztime,Zdata,ZSAChdr] = fget_sac();
    
        rx=lons;
        ry=lats;
%         times1=ZSAChdr.times.t1;
        original_sr=sr
        Opr.sr
        Ztime=((1:length(Zdata))-1)/sr;
        Zdata=decimate(Zdata,round(original_sr/Opr.sr));
        Ztime=decimate(Ztime,round(original_sr/Opr.sr));
        length(Ztime);
        Zdata=Zdata-mean(Zdata);
        
        if Opr.resample~=0
            %Ztime_re=(1:length(Ztime)*Opr.resample)/(Opr.sr*Opr.resample);
            
            Zdata=interp(Zdata,Opr.resample);
            Ztime=interp(Ztime,Opr.resample);
            
            
            
        end
        
        
        %Zdata=diff(diff(Zdata));
        if Opr.snrFilterbool==true% snr filter
            z0data=filter(B,A,Zdata);
            
            
            std(z0data((Opr.ori+Opr.snrFilter(4))*Opr.sr:(Opr.ori+Opr.snrFilter(5))*Opr.sr))/std(z0data((Opr.ori+Opr.snrFilter(6))*Opr.sr:(Opr.ori+Opr.snrFilter(7))*Opr.sr))
            if std(z0data((Opr.ori+Opr.snrFilter(4))*Opr.sr:(Opr.ori+Opr.snrFilter(5))*Opr.sr))/std(z0data((Opr.ori+Opr.snrFilter(6))*Opr.sr:(Opr.ori+Opr.snrFilter(7))*Opr.sr))>=Opr.snrFilter(1)
                continue;
                
            end
        end
        if Opr.distanceFilterbool==true% distance filter
            if distance11(ry,rx,lat0,lon0,1)*180/pi>Opr.distanceFilter(2) || distance11(ry,rx,lat0,lon0,1)*180/pi<Opr.distanceFilter(1)
                continue;
                
            end
        end
        x(count,1:length(Zdata))=Zdata;
           
        time(count,1:length(Zdata))=Ztime(1:length(Zdata));
        r(count,1)=rx;
        r(count,2)=ry;
        t1(count)=begin;
        name=stnm;
        
        nm(count,1:length(char(stnm)))=char(name);
         rdis(count)=distance11(ry,rx,lat0,lon0,1)*180/pi;
%          rx
%          ry
%          lat0
%          lon0
         [s az(count)]=vdist(ry,rx,lat0,lon0);
%         rdis(count)=ZSAChdr.evsta.dist;
%         az(count)=ZSAChdr.evsta.az;
        lat(count)=lats;
        lon(count)=lons;
        count=count+1;
    end
%     if Opr.resample~=0
%     Opr.sr=Opr.sr*Opr.resample;
%     end
    count
    nel=count-1
    B=1:nel;
    size(x)
    size(time)

ret=struct('r',r,...%the position of the stations
    'nm',char(nm),...%name of the stations
    'lat',lat,...%latitude
    'lon',lon,...%longitude
    'rdis',rdis,...%epicentral distances
    'az',az,...%azimuth of the stations
    'time',time,...% time matrix
    't1',t1,...% theoretical arrival times
    'x',x,...% original data 
    'sr',Opr.sr,...% sample rate
    'ori',Opr.ori,...% origin time
    'opr',Opr);% attached the option object to read parameters.

if Opr.nrt~=0;% n squre root of the data
    for i=1:nel
        ret.x(i,:)=sign(ret.x(i,:)).*(abs(ret.x(i,:))).^(1/Opr.nrt);
    end
end


if Opr.orderDistancebool==true% order the data according the distance
    ret=orderDistance(ret);
end
if Opr.alignbool==true% align the data
    ret=seismoAlign(ret,Opr.align);
end

if Opr.bpbool==true% bandpass ythe data
    for i=1:nel
        
        ret.x(i,:)=filter(BB,AA,ret.x(i,:));
    end
    if Opr.alignbool==true
        for i=1:nel
            
            ret.x1(i,:)=filter(BB,AA,ret.x1(i,:));
        end
    end
    
end
if Opr.dec~=0% decimate
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x01=x0;
    x=ret.x;
    x1=ret.x1;
    clear ret.x ret.x1
    for i=1:nel
          xx(i,:)=decimate(x(i,:),Opr.dec);
          xx1(i,:)=decimate(x1(i,:),Opr.dec);
    end
   ret.x=xx;
   ret.x1=xx1;
    Opr.sr=Opr.sr/Opr.dec;
    ret.sr=Opr.sr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


