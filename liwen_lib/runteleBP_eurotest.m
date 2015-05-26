% This script is used to bp 


%   libDIR = '/home/liwen/Documents/seis_research/libBP1/';
%   dataDIR0 = '/home/liwen/Documents/seis_research/raw_data/tohoku_2011_na/';
%=======what is the use of this file?==============
   %load([libDIR 'nchileCali.mat']);
%=======ret1 and ret stand for? Go to read_nchile_input_liwen=================
   %ret1=ret;
   %disp(filename);
   %load eu_testshift1.mat
   clear all;
   close all;
   load eu_testshift1.mat
   %load([dataDIR0 filename]);
%===matchshiftAll = ?===========================
   %ret1.lon=round(ret1.lon*1e4)/1e4;
   disp('test');
   ret.nf = ones(1, length(ret.lon)); % nf is the normalization factor to normalize each staions amplitude
   %[ret1 ret]=matchshiftAll(ret1,ret);
   n = ret.n;
   for i=1:n
  % ret.xori(i,:)=specshift(ret.xori(i,:),ret.sr*(ret1.recordtime(i)-ret1.recordtime(100))*24*3600);
   end
   load('ptimes.mat');
%%%%%%%%%%%%%%%%%
   ret.latrange=[-1.0, 1.0];
   ret.lonrange=[-1.0, 1.0];
   tbuff=20;
%ret.begin1=360*bh
%ret.end1=360*eh;
%bh=str2double(bh);
%eh=str2double(eh);
%if ischar(bh)
%disp(['test'])
%end
   bh = 0;
   eh = 4;
   ret.begin=max([3600*bh tbuff ]);% begin = 0 is from the first second in seismograms
   
   ret.end=min([3600*eh length(ret.xori)/ret.sr-tbuff]);
%ret.begin1
%ret.end1
   %filename;
   disp(['time begin=' num2str(ret.begin) 'end=' num2str(ret.end)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   x0 = ret.xori;
   ret = rmfield(ret,'xori');
   r = ret.r;% what is the "r" 2*n array lon lat
   lon0 = ret.lon0;
%lon0=93.063;%-70.813;
%lon0 = -70/813;
   lat0 = ret.lat0;
%lat0=2.311;%-19.642;
%lat0 = -19.642;
   sr = ret.sr;% sr = sample rate
%parr=ret.parr;% time shift at beginning here is 20s. 
   parr = 0;
   begin=ret.begin;
   over=ret.end;% time end of the time series(
%step=ret.step;% step is the time increment step, here is 1s.
step=1;
ps=40%ret.ps;% backprojection region grid row number 
qs=40%ret.qs;% backprojection region grid column number
display(['ps=' num2str(ps)]);
display(['ps=' num2str(qs)]);
uxRange=ret.latrange;% back projection region size
uyRange=ret.lonrange;% same above
fl=ret.fl;% frequency low limit
fh=ret.fh;% frequency high limit
win=ret.win;% time window in unit second 
dirname=ret.dirname;
%%%%%%%%%%%%%%%%%%%%%%%
[n m]=size(x0);
[stationNum, sampleNum]=size(x0);
display('sampleNum = ');
display( sampleNum);
ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps);
uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%saveDir=['./' dirname num2str(n) 'stations' num2str(win) 's' num2str(fl) 'HzTo' num2str(fh) 'Hz'  ] ;
%system(['mkdir ' saveDir]);
%cd(saveDir);
%fileID=fopen(['0411logfile' dataDIR '_' num2str(bh) '_' num2str(eh)],'w');
%I = findstr('.', filename);
%logname = filename(1:I(end)-1);
%fileID=fopen(['./logfile_', logname], 'w');
fileID = fopen('logfile_eutest.dat', 'w');
%fileID=fopen('logfile_eurotest.dat', 'w');
%logname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [BB,AA]=butter(4,[fl fh]/(sr/2));% butterworth filter data, the second parameter is the cutoff frequencies
   for i=1:n
       x0(i,:)=filter(BB,AA,x0(i,:));
   end
   for i=1:n 
%         x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   tlib=zeros(n,ps,qs);% what is tlib? n is the stations number.
   for p=1:ps
       for q=1:qs
           sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,tt)';
           tlib(:,p,q)=sd-mean(sd);
       end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rawData = single(x0');% return the single precision number 
   offset = repmat((0:stationNum-1)*sampleNum, [win*sr, 1]);% stationNum = n. What is win? asscociated with multitaper?  repmat is replicate and tile array
  % display(offset);
   %display(['offset size = ' offset]);
   for tl=parr+begin:step:parr+over% why add back the time shift
       th=tl+win;
       dataIndTemp = bsxfun(@plus, repmat(tl*sr, [win*sr, stationNum]), (0:win*sr-1)');% window beginning position
	   display(['dataIndTemp=' size(dataIndTemp)]);
       display(['t=' num2str(tl-parr-begin) 's']);
       Pm=zeros(ps,qs);
       Pm1=zeros(ps,qs);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	y=zeros(n,win*sr);
	y0=zeros(n,2*win*sr);
	
	for p=1:ps
		for q=1:qs
			sd=tlib(:,p,q)';
			
			dataInd = bsxfun(@plus, dataIndTemp, sd*sr);  % timeshift data ======
			 %display(dataInd);
			dataInd = floor(dataInd + offset);
            %display('dataInd = ');
          
			y = rawData(dataInd)';% Each row is time series corresponding to one station.
			% 					                    for k=1:n
			% 					                           y(k,:)=x0(k,floor((tl+sd(k))*sr):floor((tl+win+sd(k))*sr-1));
			% 					%                         y0(k,:)=specshift(x0(k,(tl-1/2*win)*sr:(tl+3/2*win)*sr-1),sd(k)*sr);
			% 					%                         y(k,:)=y0(k,1/2*win*sr:3/2*win*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
			% 					                    end
			Pm(p,q)=sum(sum(y(:,:),1).^2);% stacking    
			mati=corrcoef(y'); % mati is a 192 by 192 correlation coefficient matrix
			Pm1(p,q)=mean(mati(:)); % cross-correlation======?==============This is easy: if mati=[1,2;3,4], then mati(:)=[1,2,3,4]'
		end
	end
	tmp1=peakfit2d(real(Pm));
	bux=interp1(1:length(ux),ux,tmp1(1));
	buy=interp1(1:length(uy),uy,tmp1(2));
	maxp=interp2(1:ps,1:qs,Pm',tmp1(1),tmp1(2),'linear',0);
	maxp0=max(Pm(:));
	disp(['bm bux ' num2str(tmp1(1)) ' buy ' num2str(tmp1(2)) ' max ' num2str(maxp) ' max0 ' num2str(maxp0)]);
	
	
	%
	%
	tmp2=peakfit2d(real(Pm1));
	bux1=interp1(1:length(ux),ux,tmp2(1));
	buy1=interp1(1:length(uy),uy,tmp2(2));
	maxp1=interp2(1:ps,1:qs,Pm1',tmp2(1),tmp2(2),'linear',0);
	maxp2=interp2(1:ps,1:qs,Pm',tmp2(1),tmp2(2),'linear',0);
	disp(['cm bux ' num2str(tmp2(1)) ' buy ' num2str(tmp2(2)) ' maxp1 ' num2str(maxp1) ' maxp2 ' num2str(maxp2)]);
	
	%
	%              [maxp1 I]=max(Pm1(:));
	%              maxp2=Pm(I);
	%             fprintf(fileID,'%f %f %f %f\n',tl,bux,buy,maxp);
	fprintf(fileID,'%f %f %f %f %f %f %f %f\n',tl,bux,buy,maxp,bux1,buy1,maxp1,maxp2);
	%
	
	%         save ([  num2str(tl-parr-begin) 'smat'],'Pm');
	
end
fclose(fileID);
%save('/home/liwen/Documents/seis_research/tohoku_2011_usarray/data/parret','ret');