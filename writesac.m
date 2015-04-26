function writesac(SeisData,time,DELTA,filename,head)
% WRITESAC(SeisData,time,DELTA,filename,head)
%
% Writes time series data into SAC binary format
% 'time' contains beginning and end
% 'SeisData' contains the data stream
% 'DELTA' is the nominal increment between evenly spaced values
% May give it some header variables in a structure array
% with a selection of valid SAC variables.
% THIS TAKES PRECEDENCE OVER THE 'TIME' VARIABLE
%
% Last modified by fjsimons@alum.mit.edu, Feb. 8th, 2003
%
% Example:
%
% defval('ddir','d_20002040741');
% defval('dfil','NE_20002040741.wav');
% defval('Fs',110)
% x=readwav(fullfile('/home/fjsimons/MERMAID/',ddir,dfil));
% t=linspace(0,(length(x)-1)/Fs,length(x));
% writesac(x,t,1/110,[pref(dfil),'.sac'])
%
% Cull these values from the input or supply defaults

defval('badval',-12345);
defval('filename',[inputname(1),'.sac']);
defval('head',[])

% Get supplied header values before anything else
if ~isempty(head)
  fn=fieldnames(head);
  for index=1:size(fn,1)
    eval([fn{index} '=head.' fn{index} ';'])
  end
end

defval('NPTS',length(SeisData));
defval('time',[1 length(SeisData)])
defval('B',time(1));
defval('E',time(end));
defval('SCALE',1)
defval('INTERNAL',2)
defval('T0',badval)
defval('T1',badval)
defval('CMPAZ',0)
defval('CMPINC',0)
defval('NZYEAR',indeks(datevec(datenum(clock)),1));
defval('NZJDAY',floor(dayofyear))
defval('NZHOUR',indeks(datevec(datenum(clock)),4));
defval('NZMIN',indeks(datevec(datenum(clock)),5));
defval('NZSEC',ceil(indeks(datevec(datenum(clock)),6)));
defval('NVHDR',6)
defval('KINST','Matlab')
defval('KSTNM','Linux ')
defval('KUSER0','fjsimons')
defval('KCMPNM','single')
defval('LEVEN',1)
defval('LPSPOL',0)
defval('LCALDA',0)
defval('IFTYPE',1)
defval('DELTA',1)
defval('LOVROK',badval)
defval('IDEP',badval)
defval('IZTYPE',69)
defval('IINST',badval)
defval('ISTREG',badval)
defval('IEVREG',badval)
defval('IEVTYP',badval)
defval('IQUAL',badval)
defval('ISYNTH',badval)
defval('IMAGTYP',badval)
defval('IMAGSRC',badval)
defval('STLA',badval)
defval('STLO',badval)
defval('STEL',badval)
defval('STDP',badval)
defval('EVLA',badval)
defval('EVLO',badval)
defval('EVEL',badval)
defval('EVDP',badval)
defval('MAG',badval)

% Initialize Header
HdrFloats=repmat(badval,70,1);
HdrNhdr=repmat(badval,15,1);
HdrIhdr=repmat(badval,20,1);
HdrLhdr=repmat(badval,5,1);
HeaderStrings=repmat(sprintf('%i  ',badval),24,1);

% Assign variables to the header
HdrFloats(1)=DELTA;
HdrFloats(4)=SCALE;
HdrFloats(6)=B;
HdrFloats(7)=E;
HdrFloats(10)=INTERNAL;
HdrFloats(11)=T0;
HdrFloats(12)=T1;
HdrFloats(32)=STLA;
HdrFloats(33)=STLO;
HdrFloats(34)=STEL;
HdrFloats(35)=STDP;
HdrFloats(36)=EVLA;
HdrFloats(37)=EVLO;
HdrFloats(38)=EVEL;
HdrFloats(39)=EVDP;
HdrFloats(40)=MAG;
HdrFloats(58)=CMPAZ;
HdrFloats(59)=CMPINC;

HdrNhdr(1)=NZYEAR;
HdrNhdr(2)=NZJDAY;
HdrNhdr(3)=NZHOUR;
HdrNhdr(4)=NZMIN;
HdrNhdr(5)=NZSEC;
HdrNhdr(7)=NVHDR;
HdrNhdr(10)=NPTS;

HdrLhdr(1)=LEVEN;
HdrLhdr(2)=LPSPOL;
HdrLhdr(3)=LOVROK;
HdrLhdr(4)=LCALDA;
HdrLhdr(5)=0;

HdrIhdr(1)=IFTYPE;
HdrIhdr(2)=IDEP;
HdrIhdr(3)=IZTYPE;
HdrIhdr(5)=IINST;
HdrIhdr(6)=ISTREG;
HdrIhdr(7)=IEVREG;
HdrIhdr(8)=IEVTYP;
HdrIhdr(9)=IQUAL;
HdrIhdr(10)=ISYNTH;
HdrIhdr(11)=IMAGTYP;
HdrIhdr(12)=IMAGSRC;

HeaderStrings(24,1:length(KINST))=KINST;
HeaderStrings(1,1:length(KSTNM))=KSTNM;
HeaderStrings(18,1:length(KUSER0))=KUSER0;
HeaderStrings(21,1:length(KCMPNM))=KCMPNM;

fid=fopen(filename,'w','l');
fwrite(fid,HdrFloats,'float32');
fwrite(fid,HdrNhdr,'int32');
fwrite(fid,HdrIhdr,'int32');
fwrite(fid',HdrLhdr,'int32');
fwrite(fid,HeaderStrings','char');
fwrite(fid,SeisData,'float32');
fclose(fid);

