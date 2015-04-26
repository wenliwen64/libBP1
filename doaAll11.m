function ret=doaAll(ret,opr)
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
            'ncontour',20,...% number of contour of the diagram
            'slidingWindowViewRange',[-20 20 8],...% the seismograms plot with normalization range and ratio relative to std 
             'beamOrientation',atan2((18.788-18.6654),-(72.724-72.663))*180/pi,...% the beam direction lat/lon if you want to project to the fault line 
             'faultOrientation',atand((0.0308/0.0675)),...% the fault orientation atan(lat/lon)
             'projectionRange',[1 100],...% the range along the fault relative to the epicenter
             'saveProj',0,...% save the projection figure;
             'constantBias',-10,...%costant bias
             'saveDir','/home/lsmeng/matlab/haitiUS/gf/save/')% dir to save the figures
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
faultOrientation=opr.faultOrientation;
saveProj=opr.saveProj;
constantBias=opr.constantBias;
saveDir=opr.saveDir;
 t3=(0:length(x0)-1)/sr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(19);
hold all;
ColorSet = varycolor(20);
set(gca, 'ColorOrder', ColorSet);


kt=0;
gcfold=figure(1);
for tl=ori+begin:step:ori+over
   tl
    kt=kt+1;
    
        clear X;

        th=tl+windowLength;
        
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
            Lh=length(x0)
            y=zeros(nel,Lh);
            S=zeros(length(s),nel,nel);
            for i=1:nel
                for j=1:nel
                    %[c ph ci ]=cmtm2(x0(i,tl*sr:th*sr)/std(x0(i,tl*sr:th*sr)),x0(j,tl*sr:th*sr)/std(x0(j,tl*sr:th*sr)),2);
                    xi=(x0(i,tl*sr:th*sr-1)-mean(x0(i,tl*sr:th*sr-1)));
                    xj=(x0(j,tl*sr:th*sr-1)-mean(x0(j,tl*sr:th*sr-1)));
%                     [s c0 ph ci phi]=cmtm(xi/std(xi),xj/std(xj),1/sr,3,1);
                   [c ph ci ]=cmtm2(xi/std(xi),xj/std(xj),2);
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
                %R21b(i-fli+1)=Rxx(1,1);
                [Uv,A]=eig(Rxx);
                As=zeros(nel,nel);
                un=zeros(nel,nel);
                %us=zeros(nel,nel);
                if strcmp(Mm,'rank')
                    M=rank(Rxx);
                else
                    M=Mm;
                end
                
                
                un(:,1:nel-M)=Uv(:,1:nel-M);
                
                Un=un*un';
                
                vi=s(i);
                %wi=s1(i);
                wi=1;%ww(kkk);
                
                %t=((r(1,:)-ux(p)).^2+(r(2,:)-uy(q)).^2).^0.5/v;
                %t=t-t(1);
                for p=1:ps
                    
                    for q=1:qs
                        a=ptime(fs,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                        Pm(p,q)=Pm(p,q)+(wi*(a'*a)/(a'*Un*a));
                    end
                end
                %pause
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(method,'stacking')% stacking part
            Lh=length(x0)
            y=zeros(nel,Lh);
            for k=1:nel
                stdk(k)=std(x0(k,(tl-2*windowLength)*sr:(tl+4*windowLength)*sr));
            end
            for p=1:ps
                p
                for q=1:qs
                    
                  
                    y1=zeros(1,windowLength*sr);
                    y=zeros(nel,Lh);
                    tt=(1:length(y))/sr;
                    sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
                    sd=sd-sd(3);
                    for k=1:nel
                     
                        y(k,(tl-2*windowLength)*sr:(tl+4*windowLength)*sr)=specshift(x0(k,(tl-2*windowLength)*sr:(tl+4*windowLength)*sr),sd(k)*sr);
                        y1=y1+y(k,tl*sr:tl*sr+windowLength*sr-1)/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                        
                    end
                    Pm(p,q)=sum(y1.^2);
                    
                end
            end
        end
        stdk
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ploting the backprojection diagrams
        maxpm=max(max(Pm));
        
        
        
        for p=1:ps
            for q=1:qs
                if Pm(p,q)==maxpm
                    bux=ux(p);
                    buy=uy(q);
                end
          
                    
            end
        end
        bbux(kt)=bux;
        bbuy(kt)=buy;
        
        
        
        
        
        %scrsz = get(0,'ScreenSize');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        gcfnew=figure();
        
        
        
        
        subplot(1,2,1);
        [ch,ch]=contourf(Ym,Xm,real(Pm),ncontour);
        colorbar;
        
        
        hold on;
        plot(haiticoast(:,1),haiticoast(:,2),'white');
        plot(lon0,lat0,'r*','MarkerSize',10);
        plot(bbuy,bbux,'blackO--','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','white','MarkerSize',5);
        
        saveProj
        kt
        
        %[dist val xt yt]=tprf(Ym,Xm,Pm,buy,bux,lat0,lon0,tand(beamOrientation),faultOrientation,projectionRange);% sub rutine to comput the projection on the fault
        [dist val xt yt]=tprfa(Ym,Xm,Pm,bux,buy,lat0,lon0,beamOrientation,faultOrientation,projectionRange);
        rup(kt,1:length(val))=val/max(val);
        plot(xt,yt,'white','LineWidth',1);
        
        tt=(1:length(y))/sr;
        set(ch,'edgecolor','none');
        plot(bbux,bbuy,'white.-');
        plot(lon0-0.0675*(saveProj-1),lat0-0.0308*(saveProj-1),'r*');
        subplot(1,2,2);
        hold on;
      
        for jj=1:nel
            
            if slidingWindowViewRange(1)==slidingWindowViewRange(2)
                plot(t3,x0(jj,:)/std(x0(:))/slidingWindowViewRange(3)+(jj),'b');
            else
                plot(tt,x0(jj,:)/std(x0(jj,(ori+slidingWindowViewRange(1))*sr:(ori+slidingWindowViewRange(2))*sr))/slidingWindowViewRange(3)+(jj));
            end
        end
        
        plot([ tl tl] ,[1 nel],'r');
        plot([th th],[1 nel],'r');
        xlim([tl-5*windowLength th+7*windowLength]);
        %saveas(gcf,num2str(kt),'fig');
       % close(gcf);
        bbux(kt)=bux;
        bbuy(kt)=buy;
       % close(gcf);
        
        
        
        
        bux
        buy
        % close(gcfold);
         gcfold=gcfnew;
end
%project along the fault
%close all
gcff=figure(100);
[lx ly]=size(rup');
Xxm=zeros(lx,ly);
Yym=zeros(lx,ly);

uxxx=linspace(projectionRange(1)+constantBias,projectionRange(2)+constantBias,lx);
uyyy=linspace(0,step*ly,ly);

for q=1:ly
    Xxm(:,q)=uxxx;
end
for p=1:lx
    Yym(p,:)=uyyy';
end


pcolor(Xxm,Yym,rup');
hold on;
if saveProj ~=0
distline=distance11(lat0,lon0,lat0-0.0308*(saveProj-1),lon0-0.0675*(saveProj-1),6371);
plot(distline*[1 1],[0 100],'white','lineWidth',2);
end
if saveProj~=0
    cd(opr.saveDir);
    saveas(gcff,[num2str(saveProj) 'rupeps'],'eps');
    saveas(gcff,[num2str(saveProj) 'rupfig'],'fig');
    save([num2str(saveProj) 'rup']);
end
%shading flat;


