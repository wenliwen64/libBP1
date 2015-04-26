
clear all;
close all;
load sarret.mat;
system('ls *smat > filelist');
system('awk ''{print "cp "$1" "$1".mat"}'' filelist | sh');
figure(2);
%%%%%%%%%%%%%%%%%
sr=ret.sr;
sarr=ret.sarr;
tend=ret.end-ret.begin;
step=ret.step;
ps=ret.ps;
qs=ret.qs;
uxRange=ret.uxrange;
uyRange=ret.uxrange;
lon0=ret.lon0;
lat0=ret.lat0;
lonm=ret.lonm;
latm=ret.latm;
faultorientation=ret.faultorientation;
faultN=ret.faultN;
cct=ret.cct;
projectionRange=ret.projectionRange;
%%%%%%%%%%%%%%%%%%
ncontour=200;
t=0:step:tend-step;
Pm=zeros(ps,qs);
ux=linspace(uxRange(1),uxRange(2),ps);
uy=linspace(uyRange(1),uyRange(2),qs);
bux=zeros(length(t),10);
buy=zeros(length(t),10);
Power=zeros(length(t),10);
nw=zeros(1,length(t));

h=figure(2);
set(gcf,'Position',[100 1000 500 700]);
hold on;
load( [num2str(0) 'smat.mat']);
[c,ch]=contourf(uy,ux,real(Pm)/20,ncontour);
axis equal;
colorbar;
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);
drawnow;

f = getframe(h);
[im,map] = rgb2ind(f.cdata,256,'nodither');

clf;
imwrite(im,map,'movie.gif','DelayTime',0.50,'LoopCount',inf)
kt=0;
for j=1:length(t)
    kt=kt+1;
    
    hold on;
    title([num2str(t(j)) 's']);
    load( [num2str(t(j)) 'smat.mat']);
    
    Pm=abs(Pm);
    Pm=Pm-min(Pm(:));

          maxpm(kt)=max(max(Pm));
    [pw qw]=localMaximum1(abs(Pm));
    
    clear tmp;
    for jj=1:length(pw)
        tmp(jj)=abs(Pm(pw(jj),qw(jj)));
    end
   
    ttmp=sortrows([pw qw tmp'],-3);
    pw=ttmp(:,1);
    qw=ttmp(:,2);
    tmp=ttmp(:,3);
   nw(kt)=min([ 3 length(tmp) ]);
  
     for jj=1:nw(kt)
        if tmp(jj)<0.5*tmp(1)
            nw(kt)=jj-1;
            break;
        end
     end
     
    bux(kt,1:nw(kt))=ux(pw(1:nw(kt)));
    buy(kt,1:nw(kt))=uy(qw(1:nw(kt)));
        tmp1=peakfit2d(real(Pm));
    bux(kt,1)=interp1(1:length(ux),ux,tmp1(1));
    buy(kt,1)=interp1(1:length(uy),uy,tmp1(2));

   

%      [c,ch]=contourf(uy,ux,abs(Pm)/maxpm(kt)*10*log10(abs(Pm(pw(1),qw(1)))),ncontour);
     [c,ch]=contourf(uy,ux,abs(Pm),ncontour);


    colorbar;
    
     axis equal;
    ylim([min(ux) max(ux)]);
    xlim([min(uy) max(uy)]);
      xlabel('Se (s/km)');
      ylabel('Sn (s/km)');

    for jj=1:nw(kt)
    Power(kt,jj)=sqrt(Pm(pw(jj),qw(jj)));
    end

    set(ch,'edgecolor','none');

    chsize(20)
      drawnow;
    
    
    f=getframe(h);
    im= rgb2ind(f.cdata,map,'nodither');

    clf
    if kt==1
        imwrite(im,map,'jp_movie.gif','DelayTime',0.15,'LoopCount',inf) 
    else
   imwrite(im,map,'movie.gif','DelayTime',0.5,'WriteMode','append') 
    end
 
    dist=linspace(projectionRange(1),projectionRange(2),faultN);
    [uxx uyy val distpeak0 maxcc0]=tprfslowcali(ux,uy,Pm',bux(kt),buy(kt),lat0,lon0,latm,lonm,dist,faultorientation);
    uxxd(kt,:)=uxx;
    uyyd(kt,:)=uyy;
    rup(kt,1:length(val))=val;
    maxcc(kt)=maxcc0;
    distpeak(kt)=distpeak0;
end
save movieBP;
