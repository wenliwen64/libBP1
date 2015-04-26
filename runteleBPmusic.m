function runteleBP(DataFile)
load(DataFile);
addpath /Users/lsmeng/Dropbox/demoBP;
load ptimes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=ret.xori;
r=ret.r;
lon0=ret.lon0;
lat0=ret.lat0;
dep0=ret.dep0;
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
Nw=ret.Nw;
fs=ret.fs;
%%%%%%%%%%%%%%%%%%%%%%%

[n m]=size(x0);
ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps);
uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir=['./' dirname num2str(n) 'MUSICstations' num2str(win) 's' num2str(fl) 'HzTo' num2str(fh) 'Hz'  ] ;
system(['mkdir ' saveDir]);
cd(saveDir);
rmfield(ret,'xori');
save('parret','ret');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BB,AA]=butter(4,[fl fh]/(sr/2));
for i=1:n
    
    x0(i,:)=filter(BB,AA,x0(i,:));
end
for i=1:n 
        x0(i,:)=x0(i,:)/std(x0(i,parr*ret.sr:(parr+30)*ret.sr));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ptimesdepth=load('Ptimesdepth.mat');
rr=Ptimesdepth.dis;
dd=Ptimesdepth.dep;
tt=Ptimesdepth.Ptt;
tlib=zeros(n,ps,qs);
for p=1:ps
    
    for q=1:qs
                     
              sd=phtime2d(lat0,lon0,dep0,ux(p),uy(q),dep0,r(:,2),r(:,1),rr,dd,tt)';% along ray path
   
        tlib(:,p,q)=sd-mean(sd);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
for tl=parr+begin:step:parr+over
    
    
    th=tl+win;
    display(['t=' num2str(tl-parr-begin) 's']);
    Pm=zeros(ps,qs);
    Pw=zeros(ps,qs);

    S=cmtmall(x0(:,tl*sr:th*sr-1),Nw);
    s=linspace(0,sr-sr/((th-tl)*sr),(th-tl)*sr);
    fli=round(interp1(s,1:length(s),fl));
    fhi=round(interp1(s,1:length(s),fh));
    for i=fli:fs:fhi

        s(i)
        Pm1=zeros(ps,qs);
        Pw1=zeros(ps,qs);
        clear Uv A Un a wi;
        Rxx=zeros(n,n);
        for j=1:n
            for k=1:n
                Rxx(j,k)=S(i,j,k);
            end
        end
        [Uv,A]=eig(Rxx);
        As=zeros(n,n);
        un=zeros(n,n);
        us=zeros(n,n);
        M=rank(Rxx);
        
        
        
        un(:,1:n-M)=Uv(:,1:n-M);
        
        Un=un*un';
        for p=1:ps
            
            for q=1:qs
                
                a=exp(-1i*2*pi*s(i)*tlib(:,p,q));
                Pm1(p,q)=((a'*a)/(a'*Un*a));
                Pw1(p,q)=((a'*Rxx*a)/(a'*a));
            end
        end
        
        Pm=Pm+Pm1;
        Pw=Pw+Pw1;
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        save ([  num2str(tl-parr-begin) 'smat'],'Pm','Pw');

end
rmfield(ret,'xori');
save('parret','ret');
