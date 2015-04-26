function ret=stacking(ret)
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
           ps=20;
           qs=30;
           uxr=[-0.75 0.75];
           uyr=[-1.5 0.75]*1.5;
           lat0=18.43;
           lon0=-72.53;
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
                    
                  
                    y1=zeros(1,windowLength*sr);
                    y=zeros(n,m);
                
                    sd=phtime(1,lat0,lon0,ux(p),uy(q),lat,lon,rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:n
                       
%                         if k~=1
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(ret.x(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1=y1+y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
%                         end
                    end
                    Pm(p,q)=sum(y1.^2);
                    
                end
            end
            figure(15);
            hold on;
            
            pcolor(Ym,Xm,Pm);
            axis equal
            shading flat;
            colorbar;
            plot(haiticoast(:,1),haiticoast(:,2),'white');
            plot(lon0,lat0,'r*');
           
            ylim([ux(1) ux(end)]);
           xlim([uy(1) uy(end)]);
            hold off;