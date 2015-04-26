function runlocalBP(DataFile)
load(DataFile);
kv=@(f0,ux,uy) 2*pi*f0*[ux uy];% wave vector
v=@(r,k) exp(1i*(r*k'));%steering vector
opr.sarr=4;
x0=ret.xori;
r=ret.r;
sr=ret.sr;
sarr=ret.sarr;
begin=ret.begin;
over=ret.end;
step=ret.step;
ps=ret.ps;
qs=ret.qs;
uxRange=ret.uxrange;
uyRange=ret.uyrange;
fl=ret.fl;
fh=ret.fh;
win=ret.win;
dirname=ret.dirname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n m]=size(x0);
ux=linspace(uxRange(1),uxRange(2),ps);
uy=linspace(uyRange(1),uyRange(2),qs);
%%%%%%%%%%%%%%%%%%%%%%%
saveDir=['./' dirname num2str(n) 'stations' num2str(win) 's' num2str(fl) 'HzTo' num2str(fh) 'Hz'  ] ;
system(['mkdir ' saveDir]);
cd(saveDir);
%%%%%%%%%%%%%%%%%%%%
[BB,AA]=butter(4,[fl fh]/(ret.sr/2));
for i=1:n
    
     x0(i,:)=filter(BB,AA,x0(i,:));
end
% for i=1:n 
%         x0(i,:)=x0(i,:)/std(x0(i,sarr*ret.sr:(sarr+30)*ret.sr));
% end

for tl=sarr+begin:step:sarr+over

        th=tl+win;
        Pm=zeros(ps,qs);
        display(['t=' num2str(tl-sarr-begin) 's']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tlib=zeros(n,ps,qs);
      for p=1:ps
          
          for q=1:qs
              sd=r*[ux(p) uy(q) (1/0.3^2-ux(p)^2-uy(q)^2)^0.5]'; % assuming the near surface S wave speed of 0.3 km/s         
              tlib(:,p,q)=sd;
          end
      end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Lh=length(x0);
            y=zeros(n,win*sr);
            for p=1:ps
                for q=1:qs
                     sd=tlib(:,p,q)';

                    for k=1:n
                         y(k,:)=x0(k,floor((tl+sd(k))*sr):floor((tl+sd(k))*sr)+win*sr-1);
                    end
                    Pm(q,p)=sum(sum(y(:,:),1).^2);
                    
                end
            end
        save ([  num2str(tl-sarr-begin) 'smat'],'Pm');

end
rmfield(ret,'xori');
save('sarret','ret');
