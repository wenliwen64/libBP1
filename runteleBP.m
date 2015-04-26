function runteleBP(DataFile)
load(DataFile);
addpath ./libBP;
load ptimes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=ret.xori;
r=ret.r;
lon0=ret.lon0;
lat0=ret.lat0;
sr=ret.sr;
parr=ret.parr;
begin=ret.begin;
over=ret.end;
step=ret.step;
ps=ret.ps;
qs=ret.qs;
uxRange=ret.latrange;
uyRange=ret.lonrange;
fl=ret.fl;
fh=ret.fh;
win=ret.win;
dirname=ret.dirname;
%%%%%%%%%%%%%%%%%%%%%%%
[n m]=size(x0);
ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps);
uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir=['./' dirname num2str(n) 'stations' num2str(win) 's' num2str(fl) 'HzTo' num2str(fh) 'Hz'  ] ;
system(['mkdir ' saveDir]);
cd(saveDir);
fileID=fopen('logfile','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BB,AA]=butter(4,[fl fh]/(sr/2));
for i=1:n
    
    x0(i,:)=filter(BB,AA,x0(i,:));
end
for i=1:n 
        x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlib=zeros(n,ps,qs);
for p=1:ps
    
    for q=1:qs
                      sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,tt)';

        tlib(:,p,q)=sd-mean(sd);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
for tl=parr+begin:step:parr+over
    
 
        th=tl+win;
        display(['t=' num2str(tl-parr-begin) 's']);
        Pm=zeros(ps,qs);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  y=zeros(n,win*sr);
                  y0=zeros(n,2*win*sr);

            for p=1:ps
                
                for q=1:qs
                     sd=tlib(:,p,q)';

                    for k=1:n
                           y(k,:)=x0(k,floor((tl+sd(k))*sr):floor((tl+win+sd(k))*sr-1));
%                         y0(k,:)=specshift(x0(k,(tl-1/2*win)*sr:(tl+3/2*win)*sr-1),sd(k)*sr);
%                         y(k,:)=y0(k,1/2*win*sr:3/2*win*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                    end
                    Pm(p,q)=sum(sum(y(:,:),1).^2);
                    
                end
            end
            tmp1=peakfit2d(real(Pm));
            bux=interp1(1:length(ux),ux,tmp1(1));
            buy=interp1(1:length(uy),uy,tmp1(2));
            maxp=max(Pm(:));
            fprintf(fileID,'%f %f %f %f\n',tl,bux,buy,maxp);

        save ([  num2str(tl-parr-begin) 'smat'],'Pm');

end
fclose(fileID);
rmfield(ret,'xori');
save('parret','ret');
