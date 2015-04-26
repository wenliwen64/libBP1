function ret=doaAll1(ret,opr)
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
load worldcoast.dat
load sumatratrench.dat;
load EPGF.dat
% load /export/scratch1/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list2;
% load /export/scratch1/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list3;
% load ~/lsmeng/matlab/summatra2012/aftershock;
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
nel=length(ret.r);
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
        ux=linspace(uxRange(1)+lat0,uxRange(2)+lat0,ps);
        uy=linspace(uyRange(1)+lon0,uyRange(2)+lon0,qs);
        
        for q=1:qs
            Xm(:,q)=ux;
        end
        for p=1:ps
            Ym(p,:)=uy';
        end
        
        Pm=zeros(ps,qs);
    
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(method,'aresp')% stacking part
                Lh=length(x0);
              y=zeros(nel,Lh);
            for p=1:ps
              
                for q=1:qs
                  t=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                  
      
                  Pm(p,q)=(abs(sum(exp(-1i*2*pi*fl*t))))^2;
                   
                    
                end
            end
            Pm=Pm/max(max(Pm));
            end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if strcmp(method,'arespSum')% stacking part
                Lh=length(x0);
              y=zeros(nel,Lh);
            for p=1:ps
              
                for q=1:qs
                  t=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                  Pm(p,q)=(abs(sum(exp(-1i*2*pi*(fl+fh)*t))))^2;
                   
                    
                end
            end
            Pm=Pm/max(max(Pm));
            end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if strcmp(method,'arespDiff')% stacking part
                Lh=length(x0);
              y=zeros(nel,Lh);
            for p=1:ps
              
                for q=1:qs
                  t=phtime(fl,lat0,lon0,ux(p),uy(q),r(:,2),r(:,1),rr,ttime);
                  Pm(p,q)=(abs(sum(exp(-1i*2*pi*(fh-fl)*t))))^2;
                   
                    
                end
            end
            Pm=Pm/max(max(Pm));
            end
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
       

%         gcf5=figure(5);
            set(gcf, 'Position', [100, 100, 600, 450]); 
  set(gcf, 'PaperPositionMode', 'auto') ;
%          set(gcf,'Position',[100 100 700 700]);
         
%         
%         
%           subplot(1,2,2);
        Pm=real(Pm);
        Pm=Pm/max(Pm(:));
        [ch]=pcolor(Ym,Xm,10*log10(real(Pm)));
      
        
        colorbar;
        caxis([-20 0]);
        
        hold on;
        %%%%%%%%%%% local maximums
%           [x y]=localMaximum(real(Pm),3,true)
%           yy=interp1(1:length(uy),uy,y);
%           xx=interp1(1:length(ux),ux,x);
%           plot(yy,xx,'black.','MarkerSize',20);
        %%%%%%%%%%%%%%
%          plot(haiticoast(:,1),haiticoast(:,2),'white','lineWidth',1);
% %                   plot(worldcoast(:,1),worldcoast(:,2),'white');
% 
%          plot(EPGF(:,1),EPGF(:,2),'yellow','lineWidth',4);
%         plot(fault_trace(:,1),fault_trace(:,2),'red');
%         plot(lon0,lat0,'g*','MarkerSize',10);
% %        plot(bbuy,bbux,'blackO--','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','white','MarkerSize',5);
% load /export/scratch1/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list2;
% load /export/scratch1/lsmeng/matlab/haitiUS/paper/figures/sup/18NEICaftershock/list3;
%          scatter(list2(:,2),list2(:,1),list2(:,3)*5,'black','filled');% plot aftershocks
%          scatter(list3(:,2),list3(:,1),list3(:,3)*5,'green','filled');% plot aftershocks
        saveProj
        kt
        
%         [dist val xt yt]=tprf(Ym,Xm,Pm,buy,bux,lat0,lon0,tand(beamOrientation),faultOrientation,projectionRange);% sub rutine to comput the projection on the fault
          [dist val xt yt]=tprfa(Ym,Xm,Pm,bux,buy,lat0,lon0,beamOrientation,faultOrientation,projectionRange,N);
%                     [dist val xt yt]=tprfaline(Ym,Xm,Pm,bux,buy,lat0,lon0,beamOrientation,faultOrientation,projectionRange,N);% no projection just on the fault trace.
          rup(kt,1:length(val))=val/max(val);
          rup(:,1)
          %         plot(xt,yt,'white','LineWidth',1);
          
          tt=(1:length(y))/sr;
          set(ch,'edgecolor','none');
          %         plot(bbux,bbuy,'white.-');
          %         plot(lon0+opr.loninc*(saveProj),lat0+opr.latinc*(saveProj),'whiteo');
          plot(lon0,lat0,'whiteo');
          load worldcoast.dat
          load sumatratrench.dat;
          plot(worldcoast(:,1),worldcoast(:,2),'white');
%           scatter(aftershock(:,3),aftershock(:,2),aftershock(:,1)*4,'black','filled');
%           plot(sumatratrench(:,1),sumatratrench(:,2),'white');
          chsize(10);
%           xlabel('Longitude (^oW)');
%           ylabel('Latitude (^oN)');
          
%           axis equal;
%            title('USArray');
%                      title('VNSN');
% 
%            ylim([17.5 19.5]);
%             xlim([-73.5 -71.5]);
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
%                
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cd(opr.saveDir);
% 
%     set(gcf5, 'PaperPosition', [0 1 20 10],'PaperOrientation','landscape');
%     saveas(gcf5,[num2str(round(saveProj*10)) '-' num2str(tl) 's_snap'],'epsc');
%     saveas(gcf5,[num2str(round(saveProj*10)) '-' num2str(tl) 's_snap'],'fig');
%      print(gcf5,'-dpdf','-r600',[num2str(round(saveProj*10)) '-' num2str(tl) 's_snap']) 
%     if ismember('va',who) && ismember('vr',who)
%         save ([num2str(round(saveProj*10)) '-' num2str(tl) 's_snap_mat'],'Pm','va','vr');
%     else
      grdwrite2(uy,ux,10*log10(abs(Pm)),['arespGrdsta' num2str(nel)]);
%     end
% %         close(gcf5);
% %         
% %          set(gcf6, 'PapetoprPosition', [0 1 10 10],'PaperOrientation','landscape');
% %         saveas(gcf6,[num2str(saveProj) '-' num2str(tl) 's_snap_seis'],'epsc');
% %         saveas(gcf6,[num2str(saveProj) '-' num2str(tl) 's_snap_seis'],'fig');
% %        
% %         close(gcf6);
%         
%         bbux(kt)=bux;
%         bbuy(kt)=buy;
%        % close(gcf); sd=sd-sd(3);
%         
%         
%         
%         
%         bux;
%         buy;
% %          close(gcfold);
% %          gcfold=gcfnew;
% % valnew=val;
% end
% %project along the fault
% %close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% 
% rup=[rup(1,:); rup ;rup(end,:) ];%%%%%%%%%%add two line for visual
% 
% [lx ly]=size(rup')
% if ly~=1
%     gcff=figure(100);
%     
%     Xxm=zeros(lx,ly);
%     Yym=zeros(lx,ly);
%     
%     uxxx=-linspace(projectionRange(1)+constantBias,projectionRange(2)+constantBias,lx);
%     uyyy=linspace(begin+windowLength,over+windowLength,ly);
%     
%     for q=1:ly
%         Xxm(:,q)=uxxx;
%     end
%     for p=1:lx
%         Yym(p,:)=uyyy';
%     end
%     
%     
%     pcolor(Xxm,Yym,rup');
%     shading flat;
%     
%     
%     hold on;
% %     if saveProj ~=0
%         timeline=opr.timeline;
%         distline=distance11(lat0,lon0,lat0+opr.latinc*(saveProj),lon0+opr.loninc*(saveProj),6371);
%         distline
% %         plot(-distline*[1 1],[min(uyyy) max(uyyy)],'white','lineWidth',2);
% %          plot([min(uxxx) max(uxxx)],[1 1]*timeline,'white','lineWidth',2);
%         
% %     end
%     % if saveProj~=0
%     chsize(15);
%     pwd
%     cd(opr.saveDir);
%     saveas(gcff,[num2str(round(saveProj*10)) 'rupeps'],'epsc');
%     saveas(gcff,[num2str(round(saveProj*10)) 'rupfig'],'fig');
%     save([num2str(round(saveProj*10)) 'rup'],'rup');
%     % end
end
%  close(gcff);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ret.maxvaldist=interp1(val,dist,max(val));
% ret.rup=rup;

