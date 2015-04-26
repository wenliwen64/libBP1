% for MNM=3:15
close all;
% clearvars all -except MNM ;
clear all;
distline=zeros(1,5);
maxvaldist=zeros(1,5);
kkkt=0;
rupall=zeros(1,1);

 NN=3;

cd '/home/lsmeng/matlab/haitiUS/gf'
system(['makesynm.sh ' num2str(NN) ' ' num2str(NN*3)]);
opr=readAllSac();
opr.FileList='filelistm';
opr.lat0=18.45;
opr.lon0=-72.45;
opr.sr=40;
opr.ori=100;
opr.nrt=0;
%opr.align.timeShiftKnown=load ('/home/lsmeng/matlab/haitiUS/gf/syn/arrival000')-200;
opr.align.shiftT1bool=true;
opr.alignbool=true;
opr.bpbool=true;
opr.bp=[0.25 1];
opr.resample=0;

ret=readAllSac('/home/lsmeng/matlab/haitiUS/gf/syn/',opr);
ret0=ret;
opr.sr=ret.sr;
%ptimes=load('/home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk');
ptimes=load('ptimes');
ttime=ptimes.tt;
rr=ptimes.rr;
[nel len]=size(ret.x);
for NNN1=[1:0.2:10]
    ret=ret0;
    kkkt=kkkt+1;
   
    NNN=0;%timeshifted synthetic away from ori
    %NNN1=3;
    %close all;
t0=phtime(1,opr.lat0,opr.lon0,opr.lat0-0.01949*(NNN),opr.lon0-0.0427*(NNN),ret.r(:,2),ret.r(:,1),rr,ttime)
t01=phtime(1,opr.lat0,opr.lon0,opr.lat0-0.01949*(NNN1),opr.lon0-0.0427*(NNN1),ret.r(:,2),ret.r(:,1),rr,ttime)

%%%%%%%%%%%%%%%%%%%%%%%%%%% sin
t3=(1:length(ret.x(1,:)))/opr.sr;
% ret.x(1,:)=sin(2*pi*0.3*t3);
% ret.x1(1,:)=sin(2*pi*0.3*t3);
%%%%%%%%%%%%%%%%%%%%%% same waveformd
x2=ret.x*0;
x3=ret.x*0;
x21=ret.x*0;
x31=ret.x*0;
%tr=(NNN1-NNN)*1.875;
tr=NNN1*1.875;
sr=ret.sr;
%%%%%%%%%%%%%%%%%%%% original

% for i=1:nel
%     x2(i,:)=specshift(ret.x(i,:),-(t0(i)-t0(1))*opr.sr);
%     x21(i,:)=specshift(ret.x1(i,:),-(t0(i)-t0(1))*opr.sr);
% end
% for i=1:nel
%     x3(i,:)=specshift(ret.x(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
%     x31(i,:)=specshift(ret.x1(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
% end
% 
%     
% for i=1:nel
%    ret.x(i,:)=x2(i,:)+x3(i,:);
%    ret.x1(i,:)=x21(i,:)+x31(i,:);
% end



%%%%%%%%%%%%%%%%%%%55 reverse x2 to get x3

for i=1:nel
    x2(i,:)=specshift([zeros(1,200*sr) (ret.x(i,200*sr:270*sr)) zeros(1,length(ret.x)-270*sr-1)],-(t0(i)-t0(1))*opr.sr);
    x21(i,:)=specshift([zeros(1,200*sr) (ret.x1(i,200*sr:270*sr)) zeros(1,length(ret.x)-270*sr-1)],-(t0(i)-t0(1))*opr.sr);
end
for i=1:nel
    x3(i,:)=specshift([zeros(1,200*sr) fliplr(ret.x(i,200*sr:270*sr)) zeros(1,length(ret.x)-270*sr-1)],-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
    x31(i,:)=specshift([zeros(1,200*sr) fliplr(ret.x1(i,200*sr:270*sr)) zeros(1,length(ret.x)-270*sr-1)],-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
end

    
for i=1:nel
   ret.x(i,:)=x2(i,:)+x3(i,:);
   ret.x1(i,:)=x21(i,:)+x31(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% taper the signal so that not correlated
% for i=1:nel
%     x2(i,:)=specshift(ret.x(i,:),-(t0(i)-t0(1))*opr.sr);
%     x21(i,:)=specshift(ret.x1(i,:),-(t0(i)-t0(1))*opr.sr);
% end
% for i=1:nel
%     x3(i,:)=specshift(ret.x(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
%     x31(i,:)=specshift(ret.x1(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
% end
% 
% 
% [e v]=dpss(80*sr,2);
% e1=e(:,1)';
% e2=e(:,2)';
% 
% e11=[zeros(1,200*sr) e1 zeros(1,length(ret.x)-280*sr)];
% e22=[zeros(1,200*sr) e2 zeros(1,length(ret.x)-280*sr)];
% for i=1:nel
%    ret.x(i,:)=e11.*x2(i,:)+e22.*x3(i,:);
%    ret.x1(i,:)=e11.*x21(i,:)+e22.*x31(i,:);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%     N0=2;
%     tt0=phtime(1,opr.lat0,opr.lon0,opr.lat0-0.01949*(N0),opr.lon0-0.0427*(N0),ret.r(:,2),ret.r(:,1),rr,ttime);
%     
%     ret.opr.lat0=ret.opr.lat0-0.01949*(N0);
%     ret.opr.lon0=ret.opr.lon0-0.0427*(N0);
%     for i=1:nel
%     ret.x(i,:)=specshift(ret.x(i,:),(tt0(i))*opr.sr);
%     ret.x1(i,:)=specshift(ret.x1(i,:),(tt0(i))*opr.sr);
%     end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:nel
%     ret.x(i,:)=sign(ret.x(i,:)).*(abs(ret.x(i,:))).^(1/2);
%     ret.x1(i,:)=sign(ret.x(i,:)).*(abs(ret.x1(i,:))).^(1/2);
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%                       
% ret.x=ret.x1;
% alignopr=seismoAlign();
% alignopr.win=2;
% alignopr.ts11=199;
% alignopr.range=0.5;
% alignopr.bp=[5 10];
% ret=seismoAlign(ret,alignopr);

%%%%%%%%%%%%%%%%%
%plotSta(ret);

% plotopr=plotAll();
% plotopr.normal=[  100 160 2];
% plotopr.view=[ 100 160];
% plotopr.plotAlignbool=true;
% 
% plotAll(ret,plotopr);
%%%%%%%%%%%%%%%%%%%%%%%%
doaopr=doaAll1();
doaopr.method='music';
doaopr.freqBand=[0.3 0.3 0.3];
doaopr.ps=40;
doaopr.qs=40;
doaopr.step=tr+1;
%doaopr.begin=110+NN*3;
doaopr.begin=100;
doaopr.over=doaopr.begin+tr+1;

doaopr.windowLength=60;
doaopr.Nw=3;
%plotSpec(ret,doaopr.begin+opr.ori,doaopr.begin+opr.ori+doaopr.windowLength,1);
doaopr.slidingWindowViewRange=[100 120 4];
doaopr.uxRange=[-1.5 0.9];
doaopr.uyRange=[-1.2 1.2];
doaopr.projectionRange=[-20 100];
doaopr.beamOrientation=120;
doaopr.constantBias=0;
%doaopr.beamOrientation=-(18.5558-18.4888)/(72.4744-72.4313);
%'beamOrientation',atan2((18.4086-18.243),-(72.4084-72.3259))*180/pi,...% the beam direction lat/lon if you want to project to the fault line 
%doaopr.beamOrientation=atan2((18.4086-18.243),-(72.4084-72.3259))*180/pi;
doaopr.saveProj=NNN1;
doaopr.Mm='rank';
doaopr.saveDir=['/home/lsmeng/matlab/haitiUS/gf/save0816/reverse/'] ;

ret=doaAll1(ret,doaopr);
maxvaldist(kkkt)=ret.maxvaldist;
distline(kkkt)=sign(NNN1)*distance11(opr.lat0,opr.lon0,opr.lat0-0.01949*(NNN1),opr.lon0-0.0427*(NNN1),6371);
rupall(kkkt,1:length(ret.rup(1,:)))=ret.rup(1,:);

end
figure(101);
plot(distline,maxvaldist);
figure(102);
pcolor(rupall);
cd ..
shading flat;
     %saveas(gcf,['2source_sameGR_x' num2str(MNM)],'fig');
% end
% polyfit(distline,maxvaldist,1)
% figure(101);
% hold on;
% plot(distline,0.9585*distline+5.9039)


%%%%%%%%%%%%
% 
%                     method: 'stacking'
%                   freqBand: [0.4000 1]
%                         ps: 20
%                         qs: 20
%                    uxRange: [-3 2]
%                    uyRange: [-3 2]
%                      begin: 5
%                       over: 8
%                       step: 1
%               windowLength: 5
%                         Mm: 'rank'
%                   ncontour: 20
%     slidingWindowViewRange: [-20 20 8]
%            beamOrientation: -0.6177
%           faultOrientation: 0.3394
%            projectionRange: [1 100]