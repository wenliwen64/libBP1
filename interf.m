function ret=interf(ret)
load haiticoast.dat;
load west_hemi.dat;
load west_political.dat;
load fault_trace.dat;
load EPGF.dat
%ptimes=load('/home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk');
ptimes=load('ptimes');
ttime=ptimes.tt;
rr=ptimes.rr;      

           tl=ret.tl;
           windowLength=ret.windowLength;
           sr=ret.sr;
           [n m]=size(ret.x);
           ps=40;
           qs=60;
           uxr=[-0.75 0.75];
           uyr=[-1.5 0.75]*1.5;
           lat0=18.43;
           lon0=-72.57;
           lat=ret.lat;
           lon=ret.lon;
           
        ux=linspace(uxr(1)+lat0,uxr(2)+lat0,ps);
        uy=linspace(uyr(1)+lon0,uyr(2)+lon0,qs);
        Pm=zeros(ps,qs);
               
        Xm=zeros(ps,qs);
        Ym=zeros(ps,qs);
        Pm=zeros(ps,qs);
        
        for q=1:qs
            Xm(:,q)=ux;
        end
        for p=1:ps
            Ym(p,:)=uy';
        end
        

       
            for p=1:ps
             p
                for q=1:qs
                    
                  
                    y1=zeros(n,windowLength*sr);
                    y=zeros(n,m);
                
                    sd=phtime(1,lat0,lon0,ux(p),uy(q),lat,lon,rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:n
                       
%                         if k~=1
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(ret.x(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1(k,:)=y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
%                         end
                    end
                    mati=corrcoef(y1');
% %                      mati=corrmat(y1,10,0.5);
                    
                    Pm(p,q)=sum(sum(mati)-sum(diag(mati)));
                    
                end
            end
            figure(15);
            hold on;
%             set(gcf,'Position',[100 100 1600 600]);
%             [ch ,ch]=contourf(Ym,Xm,Pm,10);
             pcolor(Ym,Xm,Pm);
%             set(ch,'edgecolor','none');
%             pcolor(Ym,Xm,Pm);
            axis equal
            shading flat;
            colorbar;
            plot(haiticoast(:,1),haiticoast(:,2),'white');
            plot(lon0,lat0,'r*');
           
            ylim([ux(1) ux(end)]);
           xlim([uy(1) uy(end)]);
            hold off;