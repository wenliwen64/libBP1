function ret=music(ret)
load haiticoast.dat;
load west_hemi.dat;
load west_political.dat;
load fault_trace.dat;
load EPGF.dat
%ptimes=load('/home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk');
ptimes=load('ptimes');
ttime=ptimes.tt;
rr=ptimes.rr;      
Mm=ret.M;
Nw=3;
           tl=ret.tl;
           windowLength=ret.windowLength;
           sr=ret.sr;
           [nel m]=size(ret.x);
           ps=80;
           qs=80;
           uxr=ret.uxr;
           uyr=ret.uyr;
%            lat0=18.43;
%            lon0=-72.57;
           lat0=ret.lat0;
           lon0=ret.lon0;
           lat=ret.lat;
           lon=ret.lon;
           fl=ret.fl;
           fh=ret.fh;
           
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
        
        x0=ret.x;
        s=linspace(0,sr-sr/(windowLength*sr),windowLength*sr);

        fli=round(interp1(s,1:length(s),fl));
        fhi=round(interp1(s,1:length(s),fh));
        
            Lh=length(x0);
            y=zeros(nel,Lh);
            S=zeros(length(s),nel,nel);
       
            for i=1:nel
                for j=1:nel
                    %[c ph ci ]=cmtm2(x0(i,tl*sr:th*sr)/std(x0(i,tl*sr:th*sr)),x0(j,tl*sr:th*sr)/std(x0(j,tl*sr:th*sr)),2);
                    xi=(x0(i, tl*sr:tl*sr+windowLength*sr-1)-mean(x0(i,tl*sr:tl*sr+windowLength*sr-1)));
                    xj=(x0(j, tl*sr:tl*sr+windowLength*sr-1)-mean(x0(j,tl*sr:tl*sr+windowLength*sr-1)));
                    %                     [s c0 ph ci phi]=cmtm(xi/std(xi),xj/std(xj),1/sr,3,1);
                    %                     [c ph ci ]=cmtm2(xi/std(xi),xj/std(xj),Nw);
                    [c ph ci ]=cmtm2(xi,xj,Nw);
                    %                     for k=1:length(s)
                    %                         if isnan(c0(k))
                    %                             c0(k)=1;
                    %                         end
                    %
                    %                     end
                    %              c=c0.*(cosd(ph)+1i*sind(ph));
                    for k=1:length(s)
                        S(k,i,j)=c(k);
                    end
                end
            end
            
            
            
             
            for i=fli:1:fhi
                
                Pm1=zeros(ps,qs);
                % kkk=kkk+1;
                clear Uv A Un a wi;
                Rxx=zeros(nel,nel);
                %Rxx=abs(X(:,i))*abs(X(:,i))';
                %Rxx=X(:,i)*X(:,i)';
                % R21a(i-fli+1)=Rxx(1,1);
                %     Rxx=zeros(nel,nel);
                for j=1:nel
                    for k=1:nel
                        Rxx(j,k)=S(i,j,k);
                    end
                end
                %R21b(i-fli+1)=Rxx(1,1);
                [Uv,A]=eig(Rxx);
                As=zeros(nel,nel);
                un=zeros(nel,nel);
                %us=zeros(nel,nel);
                if strcmp(Mm,'rank')
                    M=rank(Rxx);
                    ['M=' num2str(M)];
                else
                    M=Mm;
                end
                
                
                 un(:,1:nel-M)=Uv(:,1:nel-M);
%                                  un(:,M+1:nel)=Uv(:,M+1:nel);
                Un=un*un';
                
                vi=s(i)
                %wi=s1(i);
                wi=1;%ww(kkk);
                
                %t=((r(1,:)-ux(p)).^2+(r(2,:)-uy(q)).^2).^0.5/v;
                %t=t-t(1);
                for p=1:ps
                    
                    for q=1:qs
                        a=ptime(vi,lat0,lon0,ux(p),uy(q),ret.r(:,2),ret.r(:,1),rr,ttime);
                       
                        Pm1(p,q)=(wi*(a'*a)/(a'*Un*a));
                    end
                end
                %pause
                  Pm1=Pm1/max(max(Pm1));
                Pm=Pm+Pm1;
%                 vi
%                 (A(nel-1,nel-1)/A(nel,nel))
                 
            end
%             Pm=real(Pm);
            figure(15);
            hold on;
%             set(gcf,'Position',[100 100 1600 600]);
%             [ch ,ch]=contourf(Ym,Xm,Pm,10);
 
             pcolor(Ym,Xm,real(Pm));
%             set(ch,'edgecolor','none');
%             pcolor(Ym,Xm,Pm);
            axis equal
            shading flat;
            colorbar;
            plot(haiticoast(:,1),haiticoast(:,2),'white');
            plot(lon0,lat0,'black*','MarkerSize',20);
           plot(west_hemi(:,1),west_hemi(:,2),'black');
           plot(west_political(:,1),west_political(:,2),'black');
            ylim([ux(1) ux(end)]);
           xlim([uy(1) uy(end)]);
            hold off;