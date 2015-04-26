function ret=doaAll2(ret,opr)
% function ret=doaAll(ret,opr)
% the backprojection function can use stacking or music

def=struct('method','stacking',...% which method to use
            'freqBand',[0.4 1 0.5],...% the frequency band to do the incoherent music
            'ps',20,...% size of grid of the backprojection diagram 
            'qs',20,...%
            'uxRange',[-3 2],...% range of the backprojection diagram in degrees
            'uyRange',[-3 2],...
            'begin',5,...%time of the first sliding window
            'over',8,...%last sliding window
            'step',1,...% step of the sliding window
            'windowLength',5,...% length of the sliding window
            'Mm','rank',...%how many eigen values to use for music signal space
            'Nw',2,...% number of taper
            'ncontour',20,...% number of contour of the diagram
            'slidingWindowViewRange',[-20 20 8],...% the seismograms plot with normalization range and ratio relative to std 
             'beamOrientation',atan2((18.788-18.6654),-(72.724-72.663))*180/pi,...% the beam direction lat/lon if you want to project to the fault line 
             'projectionRange',[1 100],...% the range along the fault relative to the epicenter
             'N',200,...%the discretization of the projection range
             'saveProj',0,...% save the projection figure;
             'timeline',0,...%draw the timeline of the second source on the color plot
             'constantBias',0,...%costant bias
             'latinc',-0.0086,...'incrementation of the faultline'
             'loninc',-0.0465,...
             'faultOrientation',atand((0.0086/0.0465)),...% the fault orientation atan(lat/lon)
             'saveDir','/home/lsmeng/matlab/haitiUS/paper/save1019/');% dir to save the figures
if ~nargin %user wants default options, give it and stop
    ret = def;
    return
elseif nargin<2
    opr=def;
end


%main settings
%%%%%%%%%%%%%%%%%%%%%%%%%%
load haiticoast.dat;
load west_hemi.dat;
load west_political.dat;
load fault_trace.dat;
load EPGF.dat
load japan_coastline.dat;
load worldcoast.dat;
% load /home/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list2;
% load /home/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list3;
%ptimes=load('/home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk');
ptimes=load('ptimes');
ttime=ptimes.tt;
rr=ptimes.rr;


sr=ret.sr;
ori=ret.ori;
x0=ret.x1;
r=ret.r;
lat0=ret.opr.lat0;
lon0=ret.opr.lon0;
[nel samlength]=size(x0);

method=opr.method;
fl=opr.freqBand(1);
fh=opr.freqBand(2);
fs=opr.freqBand(3);
ps=opr.ps;
qs=opr.qs;
uxRange=opr.uxRange;
uyRange=opr.uyRange;
begin=opr.begin;
over=opr.over;
step=opr.step;
windowLength=opr.windowLength;
Mm=opr.Mm;
ncontour=opr.ncontour;
slidingWindowViewRange=opr.slidingWindowViewRange;
beamOrientation=opr.beamOrientation;
projectionRange=opr.projectionRange;
N=opr.N;
faultOrientation=opr.faultOrientation;
saveProj=opr.saveProj;
constantBias=opr.constantBias;
saveDir=opr.saveDir;
 t3=(0:length(x0)-1)/sr;
 Nw=opr.Nw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(19);
% hold all;
% ColorSet = varycolor(20);
% set(gca, 'ColorOrder', ColorSet);


kt=0;
%gcfold=figure(1);

for tl=ori+begin:step:ori+over
   tl=round(tl);
   
    kt=kt+1;
    
        clear X;

        th=tl+windowLength;
        tl 
        th
        windowLength
         s=linspace(0,sr-sr/((th-tl)*sr),(th-tl)*sr);
               %   s=linspace(1/(th-tl),sr,(th-tl)*sr);
        
        
%         s1=linspace(0,sr,0.5*sr);
         %[s c ph ci phi]=cmtm(x0(1,tl*sr:th*sr-1),x0(2,tl*sr:th*sr-1),1/sr,3,1);
         %s=s1;
%         
%         % fspec=linspace(0,sr/2,(th-tl)*sr/2);
        fli=round(interp1(s,1:length(s),fl));
        fhi=round(interp1(s,1:length(s),fh));
       
        
        %   Rxx=zeros(nel,nel);
        %   for jj=0:9
        %       for i=1:nel
        %           X(i,:)=fft(y(i,tl*sr+jj*sr*windowLength:th*sr-1+jj*sr*windowLength));
        %
        %       end
        %       Rxx=Rxx+X(:,fli)*X(:,fli)';
        %
        %   end
        %   eig(Rxx)
        
        
       
        Xm=zeros(ps,qs);
        Ym=zeros(ps,qs);
        Pm=zeros(ps,qs);
        Pw=zeros(ps,qs);
        ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps);
        uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
        
        for q=1:qs
            Xm(:,q)=ux;
        end
        for p=1:ps
            Ym(p,:)=uy';
        end
        
        Pm=zeros(ps,qs);
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(method,'music')% music part
            Lh=length(x0);
            y=zeros(nel,Lh);
            S=zeros(length(s),nel,nel);
            for i=1:nel
                for j=1:nel
                    %[c ph ci ]=cmtm2(x0(i,tl*sr:th*sr)/std(x0(i,tl*sr:th*sr)),x0(j,tl*sr:th*sr)/std(x0(j,tl*sr:th*sr)),2);
                    xi=(x0(i,tl*sr:th*sr-1)-mean(x0(i,tl*sr:th*sr-1)));
                    xj=(x0(j,tl*sr:th*sr-1)-mean(x0(j,tl*sr:th*sr-1)));
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
            
            
            
            
            
            for i=fli:fs:fhi
                s(i)
                Pm1=zeros(ps,qs);
                 Pw1=zeros(ps,qs);
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
              
                us(:,1:M)=Uv(:,nel-M+1:nel);
%                                  un(:,M+1:nel)=Uv(:,M+1:nel);
                Un=un*un';
                sigma=0;
                for jj=1:nel-M
                sigma=sigma+1/(nel-M)*A(jj,jj);
                end
                
                
                
                 As=zeros(M,M);
                for jj=1:M
                As(jj,jj)=(A(jj,jj)-sigma);
                end
                Us=us*As*us';
%                 Us=us*(As-sigma*eye(nel))*us';
                vi=s(i);
                %wi=s1(i);
                wi=1;%ww(kkk);
                
                %t=((r(1,:)-ux(p)).^2+(r(2,:)-uy(q)).^2).^0.5/v;
                %t=t-t(1);
                for p=1:ps
                    
                    for q=1:qs
                        
                        a=ptime(vi,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                        Pm1(p,q)=(wi*(a'*a)/(a'*Un*a));
                        Pw1(p,q)=(wi*(a'*Us*a)/(a'*a));
                    end
                end
                %pause
                                b=sort(reshape(Pm1,1,ps*qs));
                                minp=mean(b(1:round(0.1*ps*qs)));
                                maxp=max(max(Pm1));
%                    Pm1=(Pm1-minp)/(maxp-minp);% normalization
                Pm=Pm+Pm1;
                Pw=Pw+Pw1;
%                 vi
%                 (A(nel-1,nel-1)/A(nel,nel))
                 
            end
            
          
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if strcmp(method,'beam')% beam part
            Lh=length(x0);
            y=zeros(nel,Lh);
            S=zeros(length(s),nel,nel);
            for i=1:nel
                for j=1:nel
                    %[c ph ci ]=cmtm2(x0(i,tl*sr:th*sr)/std(x0(i,tl*sr:th*sr)),x0(j,tl*sr:th*sr)/std(x0(j,tl*sr:th*sr)),2);
                    xi=(x0(i,tl*sr:th*sr-1)-mean(x0(i,tl*sr:th*sr-1)));
                    xj=(x0(j,tl*sr:th*sr-1)-mean(x0(j,tl*sr:th*sr-1)));
%                     [s c0 ph ci phi]=cmtm(xi/std(xi),xj/std(xj),1/sr,3,1);
                    [c ph ci ]=cmtm2(xi/std(xi),xj/std(xj),Nw);
%                       [c ph ci ]=cmtm2(xi,xj,Nw);
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
               
                
                vi=s(i);
                %wi=s1(i);
                wi=1;%ww(kkk);
       
                for p=1:ps
                    
                    for q=1:qs
                        
                        a=ptime(fs,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                        Pm(p,q)=Pm(p,q)+(wi*(a'*Rxx*a)/(a'*a));
                    end
                end
                %pause
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(method,'stacking')% stacking part
            Lh=length(x0);
            y=zeros(nel,Lh);
            for k=1:nel
                stdk(k)=std(x0(k,tl*sr:(tl+windowLength)*sr));
            end
            for p=1:ps
                p;
                for q=1:qs
                    
                  
                    y1=zeros(1,windowLength*sr);
                    y=zeros(nel,Lh);
                    tt=(1:length(y))/sr;
                    sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:nel
%                         if stdk(k)<6
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1=y1+y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
%                         end
                    end
                    Pm(p,q)=sum(abs(y1));
                    
                end
            end
            stdk
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ploting the backprojection diagrams
        if strcmp(method,'interf')% stacking part
            Lh=length(x0);
            y=zeros(nel,Lh);
            for p=1:ps
                p
                for q=1:qs
                    
                    
                    y1=zeros(nel,windowLength*sr);
                    y=zeros(nel,Lh);
                    tt=(1:length(y))/sr;
                    sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:nel
                        
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1(k,:)=y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                    end
                       mati=corrcoef(y1');
% %                      mati=corrmat(y1,10,0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% selected pairs by initial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% correlation of the arrivals
%                     if isfield(ret,'mxcr')
%                         for ii=1:nel
%                             for jj=1:nel
%                                 if ret.mxcr(ii,jj)<0.7
%                                     mati(ii,jj)=0;
%                                 end
%                             end
%                         end
%                     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Pm(p,q)=sum(sum(mati));
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % ploting the backprojection diagrams
        if strcmp(method,'cint')% stacking part
            Lh=length(x0);
            y=zeros(nel,Lh);
            for p=1:ps
                p
                for q=1:qs
                    
                    
                    y1=zeros(nel,windowLength*sr);
                    y=zeros(nel,Lh);
                    tt=(1:length(y))/sr;
                    sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:nel
                        
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1(k,:)=y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                    end
                       mati=zeros(nel,nel);
                       y11=zeros(1,nel);
                       for ii=1:nel
                           y11(ii)=sqrt(sum(y1(ii,:).^2));
                       end
% %                      mati=corrmat(y1,10,0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% selected pairs by initial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% correlation of the arrivals
                    
                    if isfield(ret,'cpair')
                        for ii=1:nel
                            for jj=1:nel
                                
                                
                                if ret.cpair(ii,jj)~=0
                                    mati(ii,jj)=sum(y1(ii,:).*y1(jj,:))/(y11(ii)*y11(jj));
                                end
                            end
                        end
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Pm(p,q)=sum(sum(mati));
                    
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if strcmp(method,'cubicstacking')% stacking part
            Lh=length(x0);
            y=zeros(nel,Lh);
            for p=1:ps
                p;
                for q=1:qs
                    
                  
                    y1=zeros(1,windowLength*sr);
                    y=zeros(nel,Lh);
                    tt=(1:length(y))/sr;
                    sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
                    sd=sd;
                    for k=1:nel
                     
                        y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                        y1=y1+nthroot(y(k,tl*sr:tl*sr+windowLength*sr-1),3);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                    end
                    Pm(p,q)=sum(abs(y1.^3));
                    
                end
            end
           end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        maxpm=max(max(Pm));
        
        
        
        for p=1:ps
            for q=1:qs
                if Pm(p,q)==maxpm
                    bux=ux(p);
                    buy=uy(q);
                    pw=p;
                    qw=q;
                end
          
                    
            end
        end
        bbux(kt)=bux;
        bbuy(kt)=buy;
        
        
        
        
        
        %scrsz = get(0,'ScreenSize');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

        gcf5=figure(5);
       
         set(gcf,'Position',[100 100 800 500]);
%         
%         
%          subplot(1,2,1);
        [ch,ch]=contourf(Ym,Xm,real(Pm),ncontour);
      
        
        colorbar;
        
        
        hold on;
        %%%%%%%%%%% local maximums
%           [x y]=localMaximum(real(Pm),3,true)
%           yy=interp1(1:length(uy),uy,y);
%           xx=interp1(1:length(ux),ux,x);
%           plot(yy,xx,'black.','MarkerSize',20);
        %%%%%%%%%%%%%%
        plot(haiticoast(:,1),haiticoast(:,2),'white');
          plot(japan_coastline(:,1),japan_coastline(:,2),'white');
          plot(worldcoast(:,1),worldcoast(:,2),'black');
        plot(EPGF(:,1),EPGF(:,2),'yellow','lineWidth',4);
        plot(fault_trace(:,1),fault_trace(:,2),'red');
        plot(lon0,lat0,'g*','MarkerSize',10);
%        plot(bbuy,bbux,'blackO--','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','white','MarkerSize',5);
% load /home/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list2;
% load /home/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list3;
%          scatter(list2(:,2),list2(:,1),list2(:,3)*5,'black','filled');% plot aftershocks
%          scatter(list3(:,2),list3(:,1),list3(:,3)*5,'green','filled');% plot aftershocks
        saveProj
        kt
        
%         [dist val xt yt]=tprf(Ym,Xm,Pm,buy,bux,lat0,lon0,tand(beamOrientation),faultOrientation,projectionRange);% sub rutine to comput the projection on the fault
          [dist val xt yt]=tprfa(Ym,Xm,Pm,bux,buy,lat0,lon0,beamOrientation,faultOrientation,projectionRange,N);
%                      [dist val xt yt]=tprfaline(Ym,Xm,Pm,bux,buy,lat0,lon0,beamOrientation,faultOrientation,projectionRange,N);% no projection just on the fault trace.
          rup(kt,1:length(val))=val/max(val);%normalized by max
        
%                     rup(kt,1:length(val))=val;
%             rup(kt,1:length(val))=(val-min(val))/(max(val)-min(val));% normalized by max and substract min
          rup(:,1)
        plot(xt,yt,'white','LineWidth',1);
        
        tt=(1:length(y))/sr;
        set(ch,'edgecolor','none');
%         plot(bbux,bbuy,'white.-');
        plot(lon0+opr.loninc*(saveProj),lat0+opr.latinc*(saveProj),'ro');
        chsize(15);
        colorbar;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  subplot(1,2,2);
%                
% %              gcf6=figure(6);
%               hold off;
%         
%                 for jj=1:nel
%         
%                     if slidingWindowViewRange(1)==slidingWindowViewRange(2)
%                        % plot(t3,x0(jj,:)/std(x0(tl*sr:th*sr))/slidingWindowViewRange(3)+(jj),'b');
%                                         plot(t3,x0(jj,:)/std(x0(1,tl*sr:th*sr))+(jj),'r');
%                     else
%                         %plot(tt,x0(jj,:)/std(x0(jj,(ori+slidingWindowViewRange(1))*sr:(ori+slidingWindowViewRange(2))*sr))/slidingWindowViewRange(3)+(jj),'g');
%                         plot(t3,x0(jj,:)/std(x0(1,(tl)*sr:(th)*sr))/slidingWindowViewRange(3)+(jj),'b');
%                     end
% %                     text(tl,jj,num2str(stdk(jj)));
%                     hold on;
%                 end
%                 title(['start at' num2str(tl) 's']);
%                 plot([ tl tl] ,[1 nel],'r');
%                 plot([th th],[1 nel],'r');
%                 xlim([tl-5 th+5]);
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(opr.saveDir);

    set(gcf5, 'PaperPosition', [0 1 20 10],'PaperOrientation','landscape');
        saveas(gcf5,[  num2str(tl) 's_snap'],'epsc');
        saveas(gcf5,[ num2str(tl) 's_snap'],'fig');
             if strcmp(method,'music')
                 pwm=Pw(pw,qw);
        save ([  num2str(tl) 's_snap_mat'],'Pm','fli','fhi','pw','qw','Pw','pwm','Rxx','ux','uy');
             else
        save ([  num2str(tl) 's_snap_mat'],'Pm');
             end
        close(gcf5);
%         
%          set(gcf6, 'PapetoprPosition', [0 1 10 10],'PaperOrientation','landscape');
%         saveas(gcf6,[num2str(saveProj) '-' num2str(tl) 's_snap_seis'],'epsc');
%         saveas(gcf6,[num2str(saveProj) '-' num2str(tl) 's_snap_seis'],'fig');
%        
%         close(gcf6);
        
        bbux(kt)=bux;
        bbuy(kt)=buy;
       % close(gcf); sd=sd-sd(3);
        
        
        
        
        bux;
        buy;
%          close(gcfold);
%          gcfold=gcfnew;
% valnew=val;
end
%project along the fault
%close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

rup=[rup(1,:); rup ;rup(end,:) ];%%%%%%%%%%add two line for visual
  ret.rup=rup;
[lx ly]=size(rup')
if ly~=1
    gcff=figure(100);
    
    Xxm=zeros(lx,ly);
    Yym=zeros(lx,ly);
    
    uxxx=linspace(projectionRange(1)+constantBias,projectionRange(2)+constantBias,lx);
    uyyy=linspace(begin,over,ly);
    
    for q=1:ly
        Xxm(:,q)=uxxx;
    end
    for p=1:lx
        Yym(p,:)=uyyy';
    end
    
    
    pcolor(Xxm,Yym,rup');
    shading flat;
    
    
    hold on;
%     if saveProj ~=0
        timeline=opr.timeline;
        distline=distance11(lat0,lon0,lat0+opr.latinc*(saveProj),lon0+opr.loninc*(saveProj),6371);
        distline
%         plot(-distline*[1 1],[min(uyyy)
%         max(uyyy)],'white','lineWidth',2);
%          plot([min(uxxx) max(uxxx)],[1 1]*timeline,'white','lineWidth',2);
        
%     end
    % if saveProj~=0
    chsize(15);
    pwd
    cd(opr.saveDir);
    saveas(gcff, 'rupeps','epsc');
    saveas(gcff,'rupfig','fig');
%     save(rup,'rup');
    % end
    caxis([0 1]);
end
 close(gcff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ret.maxvaldist=interp1(val,dist,max(val));
ret.rup=rup;

