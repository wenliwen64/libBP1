function runteleBPcontHoff1(dataDIR)
%cd(['/u/scratch/l/lsmeng/cont/' dataDIR]);
load('/u/scratch/l/lsmeng/cont/libBP1/nchileCali.mat');
ret1=ret;
disp(dataDIR);
load(['/u/scratch/l/lsmeng/cont/' dataDIR '/nchileinput.mat']);
ret.nf=ones(1,length(ret.lon));
[ret1 ret]=matchshiftAll(ret1,ret);
n=length(ret.lon);
for i=1:n
    ret.xori(i,:)=specshift(ret.xori(i,:),ret.sr*(ret1.recordtime(i)-ret1.recordtime(100))*24*3600);
end
load('/u/scratch/l/lsmeng/cont/libBP1/ptimes.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=ret.xori;
ret=rmfield(ret,'xori');
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
%saveDir=['./' dirname num2str(n) 'stations' num2str(win) 's' num2str(fl) 'HzTo' num2str(fh) 'Hz'  ] ;
%system(['mkdir ' saveDir]);
%cd(saveDir);
fileID=fopen(['logfile' dataDIR],'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BB,AA]=butter(4,[fl fh]/(sr/2));
for i=1:n
    
    x0(i,:)=filter(BB,AA,x0(i,:));
end
for i=1:n 
%         x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
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
        Pm=zeros(ps,qs);      
        Pm1=zeros(ps,qs);      
for tl=parr+begin:step:parr+over
    
 
        th=tl+win;
        display(['t=' num2str(tl-parr-begin) 's']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  y=zeros(n,win*sr);
%                   y0=zeros(n,2*win*sr);

            for p=1:ps
                
                for q=1:qs
                     sd=tlib(:,p,q)';

                    for k=1:n
                           y(k,:)=x0(k,floor((tl+sd(k))*sr):floor((tl+win+sd(k))*sr-1));
                           
%                         y0(k,:)=specshift(x0(k,(tl-1/2*win)*sr:(tl+3/2*win)*sr-1),sd(k)*sr);
%                         y(k,:)=y0(k,1/2*win*sr:3/2*win*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                    end
                    Pm(p,q)=sum(sum(y(:,:),1).^2);
%                        mati=corrcoef(y');
%                     Pm1(p,q)=sum(sum(mati));
                end
            end
            tmp1=peakfit2d(real(Pm));
            bux=interp1(1:length(ux),ux,tmp1(1));
            buy=interp1(1:length(uy),uy,tmp1(2));
            maxp=interp2(1:ps,1:qs,Pm',tmp1(1),tmp1(2),'linear',0);
            %maxp1=max(Pm(:));
            disp(['bm bux ' num2str(tmp1(1)) ' buy ' num2str(tmp1(2)) ' max ' num2str(maxp)]);
            

% 
%             
%             tmp2=peakfit2d(real(Pm1));
%             bux1=interp1(1:length(ux),ux,tmp2(1));
%             buy1=interp1(1:length(uy),uy,tmp2(2));
%                        disp(['bm bux' num2str(tmp2(1)) 'buy' num2str(tmp2(2))]);
% %            maxp1=interp2(1:ps,1:qs,Pm1,tmp2(1),tmp2(2),'linear',0);
% %            maxp2=interp2(1:ps,1:qs,Pm,tmp2(1),tmp2(2),'linear',0);
%             
%              [maxp1 I]=max(Pm1(:));
%              maxp2=Pm(I);
             fprintf(fileID,'%f %f %f %f\n',tl,bux,buy,maxp);
% %             fprintf(fileID,'%f %f %f %f %f %f %f %f\n',tl,bux,buy,maxp,bux1,buy1,maxp1,maxp2);
% 

%         save ([  num2str(tl-parr-begin) 'smat'],'Pm');

end
fclose(fileID);
rmfield(ret,'xori');
save('parret','ret');
