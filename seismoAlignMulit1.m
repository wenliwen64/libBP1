  function ret = seismoAlignMulit(ret,align)
% function ret = seismoAlign(ret,align)
% function align the data matrix ret.x with first arrivals by
% crosscorrelate all the seismogram with a reference seismogram

def=struct('win',9,...% length of the window for the cross correlation
    'lt',100,...% range*2 is devided by lt as the step of the cross correlation
    'ts11',140,...% starting time of the window
    'range',10,...% plus and minus range of the possible time shift
    'cutoff',0,...% a correlation coefficient cut off below which will be throw away
    'refst',1,...% the number of the reference seismograms
    'bpbool',0,...%whether to do the cross correlation within a frequency band
    'bp',[0.5 2],...% the frequency band
     'timeShiftKnown',0,...,%input know timeshift default 0:unknown timeshift. a vector known timeshift suppress the timeshift caculated.
     'shiftT1bool',false);% use theoretical arrival time to align data
                         % can be get by gawk '{print "saclst t1 f 1_*"$1}' filelist1 | sh > arrival;
% if ~nargin %user wants default options, give it and stop
%     ret= def;
%     return
% elseif nargin<2
%     align=def;
% end

%%%%%%%%%%%% load data when testing
% load ~/matlab/haitiUS/savetmp/retmain;
% align=alignopr;
%%%%%%%%
[nel lengthSam]=size(ret.x);
x=ret.x;
x0=ret.x;
sr=ret.sr;
[BB,AA]=butter(4,align.bp/(sr/2));

if align.bpbool==true
    for i=1:nel
        x(i,:)=filter(BB,AA,x(i,:));
    end
   
end

%parameter settings

win=align.win;%9
lt=align.lt;%200
ts11=align.ts11;%ori
range=align.range;%10
cutoff=align.cutoff;%0.85
refst=align.refst;%1
  timeshift=zeros(1,nel);
   t3=(0:length(x)-1)/sr;
 timeShiftKnown=align.timeShiftKnown;
 shiftT1bool=align.shiftT1bool;
 B=1:nel;
 
 if timeShiftKnown(1)~=0% unknowtimeshift need to caculate
     
     timeshift=timeShiftKnown;
 elseif shiftT1bool==true;
     timeshift=ret.t1-200;
 else
    
     
     tou=-range:1/sr:range;%linspace(-range,range,2*range*sr);% possible time shift range
     lt=length(tou);
     xcr=zeros(nel,nel,lt);% all correlation coeff
     mxcr=zeros(nel,nel);% maxium cross corelation coeff
     timeshift=zeros(nel,nel);
     for k=1:nel
         
         xx0=x(k,ts11*sr:(ts11+win)*sr);%reference seismograms
         xx02=sum(xx0.^2);
         for i=1:nel
             k
             i
%              ['k= ' num2str(k) ' i= ' mum2str(i)]
             for j=1:lt
                 clear tip xx;
                 ts22=ts11+tou(j);
%                  tip=linspace(ts22,ts22+win,win*sr+1);
                 xx=x(i,ts22*sr:(ts22+win)*sr);
                 a(j)=sum(xx0.*xx)/sqrt(xx02*sum(xx.^2));
                 xcr(k,i,j)=sum(xx0.*xx)/sqrt(xx02*sum(xx.^2));
%               1                     xcr(k,i,j)=-sum(abs(xx0-xx));

             end
           
         end
         
         
         
         
         for i=1:nel
            [ mma, maxi]=max(xcr(k,i,:));
            if maxi>1 && maxi<lt
            xc1 = 0.5*( xcr(k,i,maxi+1)-xcr(k,i,maxi-1) );
            xc2 = xcr(k,i,maxi-1)-2*xcr(k,i,maxi)+xcr(k,i,maxi+1);
            maxi = maxi-xc1/xc2;
            mma = mma-0.5*xc1*xc1/xc2;
            timeshift(k,i)=interp1(1:lt,tou,maxi);
            mxcr(k,i)=mma;
            end
         end
     
     end
    
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%make the matrix symmetric
 for i=1:nel
     for j=1:nel
         if i>j
             if mxcr(i,j)>mxcr(j,i);
                 mxcr(j,i)=mxcr(i,j);
                 timeshift(j,i)=-timeshift(i,j);
             else
                 mxcr(i,j)=mxcr(j,i);
                 timeshift(i,j)=-timeshift(j,i);
             end
         end
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X=zeros(nel*(nel-1)/2+1,nel);
X=zeros(1,nel);
Y=zeros(1,1);
W=zeros(1,1);
count=1;
for i=1:nel
    for j=1:nel
        
         if i~=j
            W(count)=1/(1-mxcr(i,j));
            Y(count)=W(count)*timeshift(i,j);
            X(count,i)=-1*W(count);
            X(count,j)=1*W(count);
            
            count=count+1;
         end
    end
end

X(count,:)=ones(1,nel);
Y(count)=0;
[timeshift1 stats]=robustfit(X,Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     shifting the data with y2;
    x1=zeros(nel,length(x0));
    for i=1:nel
          x1(i,:)=specshift(x0(i,:),timeshift1(i+1)*sr);
    end
    %%%%%%%%%%%%%%%%%%%%%
    if isfield(ret,'xori')
        for i=1:nel
            ret.xori(i,:)=specshift(ret.xori(i,:),timeshift1(i+1)*sr);
            %                      x1(i,:)=mxcr(i)/abs(mxcr(i))*specshift(x0(i,:),timeshift(i)*sr);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B=find(B(:)~=0)
%     ret.x=x0(B,:);
    ret.r=ret.r(B,:);
    ret.nm=ret.nm(B,:);
    ret.rdis=ret.rdis(B);
%     ret.az=ret.az(B);
    ret.x1=x1(B,:);
     ret.x=x1(B,:);
%     ret.t1=ret.t1(B);% the aligned the seismograms
    ret.timeshift=timeshift;% the timeshift 
    ret.mxcr=mxcr;
    ret.timeshift1=timeshift1;
    ret.mxcrstats=stats;
    %%%%%%%%%%%%%%%%%%%%%%
    
    figure(200);
    ret.order=corrmapmat(mxcr,ret.nm);