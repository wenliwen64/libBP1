function ret=intf(ret)
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
        
        P=zeros(n,windowLength*sr);
       for i=1:n
           P(i,:)=fft(ret.x(i,tl*sr:tl*sr+windowLength*sr-1));
       end
       F=linspace(0,sr,windowLength*sr);
       fl=ret.fl;
       fh=ret.fh;
       fli=round(interp1(F,1:length(F),fl));
       fhi=round(interp1(F,1:length(F),fh));
       xd=ret.xd;
       distx=ones(length(fli:fhi),n,n);
       
       for i=fli:fhi
                        for ii=1:n
                            for jj=1:n
                                if distance22(ret.lat(ii),ret.lon(ii),ret.lat(jj),ret.lon(jj)) >xd/F(i)
                                distx(i-fli+1,ii,jj)=0;
                                end
                            end
                        end
       end
            for p=1:ps
             p
                for q=1:qs
                    
                  
                 
                    sd=phtime(1,lat0,lon0,ux(p),uy(q),lat,lon,rr,ttime)';
                    ct=0;
                    for i=fli:fhi
                        for ii=1:n
                            for jj=1:ii
                                if distx(i-fli+1,ii,jj) ==1
                                ct=ct+P(ii,i)*P(jj,i)'*exp(-1i*2*pi*F(i)*(sd(jj)-sd(ii)));
                                end
                            end
                        end
                    end
                  

                    
                    Pm(p,q)=abs(ct);
                    
                end
            end
            ret.Pm=Pm;
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