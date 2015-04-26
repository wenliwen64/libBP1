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



%%%%%%%%%%%%%%%%%%%%%%%%%%%main settings

load worldcoast.dat;

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
saveDir=opr.saveDir;
t3=(0:length(x0)-1)/sr;
Nw=opr.Nw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



kt=0;


for tl=ori+begin:step:ori+over % looping time
    tl=round(tl);
    
    kt=kt+1;
    
    clear X;
    
    th=tl+windowLength;
    tl
    th
    windowLength
    s=linspace(0,sr-sr/((th-tl)*sr),(th-tl)*sr);

    fli=round(interp1(s,1:length(s),fl));
    fhi=round(interp1(s,1:length(s),fh));
    

    
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% computing travel times
    tlib=zeros(nel,ps,qs);
    for p=1:ps
        
        for q=1:qs
            sd=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime)';
            tlib(:,p,q)=sd;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method,'music')% Multitaper-Music back projection
        Lh=length(x0);
        y=zeros(nel,Lh);
        
        S=cmtmall(x0(:,tl*sr:th*sr-1),Nw);
        
        
        
        
        
        for i=fli:fs:fhi
            s(i)
            Pm1=zeros(ps,qs);
            Pw1=zeros(ps,qs);
            clear Uv A Un a wi;
            Rxx=zeros(nel,nel);
   
            for j=1:nel
                for k=1:nel
                    Rxx(j,k)=S(i,j,k);
                end
            end
            [Uv,A]=eig(Rxx);
            As=zeros(nel,nel);
            un=zeros(nel,nel);
            if strcmp(Mm,'rank')
                M=rank(Rxx);
                ['M=' num2str(M)];
            else
                M=Mm;
            end
            
            
            un(:,1:nel-M)=Uv(:,1:nel-M);
            us(:,1:M)=Uv(:,nel-M+1:nel);
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
            vi=s(i);
            wi=1;
            for p=1:ps
                
                for q=1:qs
                    
                    a=exp(-1i*2*pi*vi*tlib(:,p,q));
                    Pm1(p,q)=(wi*(a'*a)/(a'*Un*a));% MUSIC spectrum
                    Pw1(p,q)=(wi*(a'*Rxx*a)/(a'*a)); % beamforming power
                end
            end
            b=sort(reshape(Pm1,1,ps*qs));
            minp=mean(b(1:round(0.1*ps*qs)));
            maxp=max(max(Pm1));%old normalization
            %                     Pm1=(Pm1-minp)/(maxp-minp);% normalization
            Pm1=Pm1/max(max(Pm1));
            
            Pm=Pm+Pm1;
            Pw=Pw+Pw1;

            
        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method,'beam')% beam part
        Lh=length(x0);
        y=zeros(nel,Lh);
        S=zeros(length(s),nel,nel);
        
        S=cmtmall(x0(:,tl*sr:th*sr-1),Nw);
        
        
        
        
        
        for i=fli:1:fhi
            clear Uv A Un a wi;
            Rxx=zeros(nel,nel);
            
            for j=1:nel
                for k=1:nel
                    Rxx(j,k)=S(i,j,k);
                end
            end
            
            
            vi=s(i);
            wi=1;
            
            for p=1:ps
                
                for q=1:qs
                    
                    a=exp(-1i*2*pi*vi*tlib(:,p,q));
                    Pm(p,q)=Pm(p,q)+(wi*(a'*Rxx*a)/(a'*a));
                end
            end
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method,'stacking')% stacking part
        Lh=length(x0);
        y=zeros(nel,Lh);
        for p=1:ps
            p
            for q=1:qs
                
                
                y1=zeros(1,windowLength*sr);
                y=zeros(nel,Lh);
                tt=(1:length(y))/sr;
                sd=tlib(:,p,q)';
                if isfield(ret,'refstation')
                    if ret.refstation==1
                        sd=sd-sd(1);
                    elseif ret.refstation==2
                        sd=sd-mean(sd);
                        
                    end
                end
                for k=1:nel
                    
                    y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                    y1=y1+y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                end
                Pm(p,q)=sum(y1.^2);
                
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(method,'interf')% stacking part
        Lh=length(x0);
        y=zeros(nel,Lh);
        for p=1:ps
            p
            for q=1:qs
                
                
                y1=zeros(nel,windowLength*sr);
                y=zeros(nel,Lh);
                tt=(1:length(y))/sr;
                sd=tlib(:,p,q)';
                %                     sd=sd-sd(3);
                for k=1:nel
                    
                    y(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr)=specshift(x0(k,(tl-1*windowLength)*sr:(tl+1*windowLength)*sr),sd(k)*sr);
                    y1(k,:)=y(k,tl*sr:tl*sr+windowLength*sr-1);%/std(y(k,tl*sr:tl*sr+windowLength*sr-1));
                end
                mati=corrcoef(y1');
                
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
                sd=tlib(:,p,q)';
                if isfield(ret,'refstation')
                    if ret.refstation==1
                        sd=sd-sd(1);
                    elseif ret.refstation==2
                        sd=sd-mean(sd);
                        
                    end
                end
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
    

    
    gcf5=figure(5);
    
    set(gcf,'Position',[100 100 800 500]);
 
    [ch,ch]=contourf(Ym,Xm,real(Pm),ncontour);
    
    
    colorbar;
    
    
    hold on;


    plot(worldcoast(:,1),worldcoast(:,2),'black');

    plot(lon0,lat0,'g*','MarkerSize',10);


 

    
    tt=(1:length(y))/sr;
    set(ch,'edgecolor','none');
    chsize(15);
    xlim([min(uy) max(uy)]);
    ylim([min(ux) max(ux)]);
    colorbar;
    
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
    bbux(kt)=bux;
    bbuy(kt)=buy;
    
end

