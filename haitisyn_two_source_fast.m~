% for MNM=3:15
close all;
% clearvars all -except MNM ;
clear all;
distline=zeros(1,5);
maxvaldist=zeros(1,5);
kkkt=0;
rupall=zeros(1,1);
%timeshifted synthetic away from ori
%NNN1=3;
%close all;
cd '/home/lsmeng/matlab/haitiUS/gf'
 NN=3;
system(['makesynm.sh ' num2str(NN) ' ' num2str(NN*3)]);


opr=readAllSac();
opr.FileList='filelistm';
opr.lat0=18.45;
opr.lon0=-72.45;
opr.latinc=-0.0086;
opr.loninc=-0.0465;
opr.sr=40;
opr.ori=100;
opr.nrt=0;
%opr.align.timeShiftKnown=load ('/home/lsmeng/matlab/haitiUS/gf/syn/arrival000')-200;
opr.align.shiftT1bool=true;
opr.alignbool=true;
opr.bpbool=true;
opr.bp=[0.2 1];
opr.resample=0;

ret1=readAllSac('/home/lsmeng/matlab/haitiUS/gf/syn/',opr);
ret=ret1;
opr.sr=ret.sr;
%ptimes=load('/home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk');
ptimes=load('ptimes');
ttime=ptimes.tt;
rr=ptimes.rr;
[nel len]=size(ret.x);


cd '/home/lsmeng/matlab/haitiUS/gf'
 NN=3;
system(['makesynm1.sh ' num2str(NN) ' ' num2str(NN*3)]);

ret2=readAllSac('/home/lsmeng/matlab/haitiUS/gf/syn1/',opr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% randomlize the phase of the second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b0=200;
e0=270;
figure(7);
ttt=(1:length(ret.x))/ret.sr;
% plot(ttt,ret1.x(1,:),ttt,ret2.x(1,:),'r');
tmp=zeros(1,1);
len=length(ret.x1(1,b0*ret.sr:e0*ret.sr));
r=2*pi*rand(1,(len-1)/2);
for i=1:nel
    tmp(i,1:len)=ret2.x1(i,b0*ret.sr:e0*ret.sr);
    tmpf=fft(tmp(i,:));
    tmpa=abs(tmpf);
    tmpph=tmpf./tmpa;
    
  r=2*pi*rand(1,(len-1)/2);% making noise
    
    tmpph(2:(length(tmpph)-1)/2+1)=cos(r)+1i*sin(r);
    tmpph((length(tmpph)-1)/2+2:length(tmpph))=fliplr(tmpph(2:(length(tmpph)-1)/2+1))';
    tmp(i,1:len)=ifft(tmpph.*tmpa);
    ret2.x1(i,b0*ret.sr:e0*ret.sr)=tmp(i,1:len);
end
ttt1=ttt(b0*ret.sr:e0*ret.sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% use the taper as the signal
% 
%     N=256*8;
% NW=30/2;
% [E,V] = dpss(N,NW);
% sr=40;
% t=(1:N)/sr;
% figure(20);
% hold on;
% for i=NW*2-2:NW*2
%     plot(t,E(:,i),'color',rand(1,3));
% end
% 
% figure(21);
%  hold on;
% for i=NW*2-2:NW*2
% [Pxx,w] = pmtm(E(:,i),[]);
% plot(w/(2*pi)*sr,Pxx,'color',rand(1,3));
%  xlim([0 2])
% end
% for i=1:nel
%     %         ret1.x1(i,:)=filter(BB,AA,ret1.x1(i,:));
%     %         ret2.x1(i,:)=filter(BB,AA,ret2.x1(i,:));
%     ret2.x1(i,b0*ret.sr:b0*ret.sr+N-1)=E(:,NW*2-1)'*1e6;
%     ret1.x1(i,b0*ret.sr:b0*ret.sr+N-1)=E(:,NW*2)'*1e6;
% end
%     



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% filter to narrow band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% or use harmonics as
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the signal
%  [BB,AA]=butter(4,[0.25 0.35]/(ret.sr/2)); 
%     for i=1:nel
% %         ret1.x1(i,:)=filter(BB,AA,ret1.x1(i,:));
% %         ret2.x1(i,:)=filter(BB,AA,ret2.x1(i,:));
%         ret2.x1(i,:)=sin(2*pi*0.002*ttt+1).*sin(2*pi*0.3*ttt);
%         ret1.x1(i,:)=cos(2*pi*0.003*ttt+1).*sin(2*pi*0.3*ttt);
%     end
      
    
    
    
    % ttt=(1:length(ret.x))/ret.sr;
    
    
    
    
    
    
    
    
    
    


 ttt=(1:length(ret.x))/ret.sr;
figure(7);
plot(ttt,ret1.x1(1,:),ttt,ret2.x1(1,:),'r');
% xlim([b0 e0]);
figure(8);






NNNb=5;
NNNe=6;
for NNN1=[NNNb:0.1:NNNe]
    
    kkkt=kkkt+1;
   
    NNN=0;
    t0=phtime(1,opr.lat0,opr.lon0,opr.lat0+opr.latinc*(NNN),opr.lon0+opr.loninc*(NNN),ret1.r(:,2),ret1.r(:,1),rr,ttime);
    t01=phtime(1,opr.lat0,opr.lon0,opr.lat0+opr.latinc*(NNN1),opr.lon0+opr.loninc*(NNN1),ret1.r(:,2),ret1.r(:,1),rr,ttime);
    %%%%%%%%%%%%%%%%%%%%%%%%%%% sin
    t3=(1:length(ret1.x(1,:)))/opr.sr;
    % ret.x(1,:)=sin(2*pi*0.3*t3);
    % ret.x1(1,:)=sin(2*pi*0.3*t3);
    %%%%%%%%%%%%%%%%%%%%%% same waveformd
    x2=ret.x*0;
    x3=ret.x*0;
    x21=ret.x*0;
    x31=ret.x*0;
    %tr=(NNN1-NNN)*1.875;
     tr=NNN1*1.875;
    %tr=0;
    for i=1:nel
        x2(i,:)=specshift(ret1.x(i,:),-(t0(i)-t0(1))*opr.sr);
        x21(i,:)=specshift(ret1.x1(i,:),-(t0(i)-t0(1))*opr.sr);
    end
    for i=1:nel
        x3(i,:)=specshift(ret1.x(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
        x31(i,:)=specshift(ret1.x1(i,:),-(t01(i)-t01(1))*opr.sr-tr*opr.sr);
    end
    
    
    % [e v]=dpss(80*sr,2);
    % e1=e(:,1);
    % e2=e(:,2);
    %
    % e11=[zeros(1,200*sr) e1 zeros(1,length(ret.x)-200*sr)]
    for i=1:nel
        ret.x(i,:)=x2(i,:)+x3(i,:);
        ret.x1(i,:)=x21(i,:)+x31(i,:);
%         ret.x(i,:)=x3(i,:);
%         ret.x1(i,:)=x31(i,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %     N0=2;
    %     tt0=phtime(1,opr.lat0,opr.lon0,opr.lat0+opr.latinc*(N0),opr.lon0+opr.loninc*(N0),ret.r(:,2),ret.r(:,1),rr,ttime);
    %
    %     ret.opr.lat0=ret.opr.lat0+opr.latinc*(N0);
    %     ret.opr.lon0=ret.opr.lon0+opr.loninc*(N0);
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
%     plotSta(ret);
    
    plotopr=plotAll();
    plotopr.normal=[  100 160 2];
    plotopr.view=[ 100 160];
    plotopr.plotAlignbool=true;
    
%     plotAll(ret,plotopr);
    
    doaopr=doaAll1();
    doaopr.method='music';
    doaopr.freqBand=[0.3 0.3 0.3];
    doaopr.ps=40;
    doaopr.qs=40;
    doaopr.step=tr/3;
    %doaopr.begin=110+NN*3;
    doaopr.begin=100;
    doaopr.over=doaopr.begin+tr;
    
    doaopr.windowLength=60;
    doaopr.Nw=3;
%     plotSpec(ret,doaopr.begin+opr.ori,doaopr.begin+opr.ori+doaopr.windowLength,1);
    doaopr.slidingWindowViewRange=[100 120 4];
    doaopr.uxRange=[-1.5 0.9];
    doaopr.uyRange=[-1.2 1.2];
    doaopr.projectionRange=[-20 100];
    doaopr.beamOrientation=123;
    doaopr.constantBias=0;
    doaopr.latinc=opr.latinc;
    doaopr.loninc=opr.loninc;
    
    
    %doaopr.beamOrientation=-(18.5558-18.4888)/(72.4744-72.4313);
    %'beamOrientation',atan2((18.4086-18.243),-(72.4084-72.3259))*180/pi,...% the beam direction lat/lon if you want to project to the fault line
    %doaopr.beamOrientation=atan2((18.4086-18.243),-(72.4084-72.3259))*180/pi;
    doaopr.saveProj=NNN1;
    doaopr.Mm='rank';
    doaopr.saveDir=['/home/lsmeng/matlab/haitiUS/paper/distanceVSsepration_music1to1/'] ;
    
    ret=doaAll1(ret,doaopr);
    NNN1
    maxvaldist(kkkt)=ret.maxvaldist;
    distline(kkkt)=sign(NNN1)*distance11(opr.lat0,opr.lon0,opr.lat0+opr.latinc*(NNN1),opr.lon0+opr.loninc*(NNN1),6371);
    %%%%%%%%%%%%%%%%%%%
    [lens lenrup]=size(ret.rup);
    rupmax=zeros(1,lenrup);
    for j=1:lenrup
        rupmax(j)=max(ret.rup(:,j));
    end
%     rupall(kkkt,1:length(ret.rup(1,:)))=ret.rup(1,:);
    rupall(kkkt,1:length(ret.rup(1,:)))=rupmax;
    %%%%%%%%%%%%%%%%%%%%%%%%
end
figure(101);
plot(distline,maxvaldist);
gcf102=figure(102);

[lx ly]=size(rupall);
Xxm=zeros(lx,ly);
Yym=zeros(lx,ly);

uxxx=linspace(NNNb,NNNe,lx)*5;
uyyy=linspace(-20,100,ly);

for q=1:ly
    Xxm(:,q)=uxxx;
end
for p=1:lx
    Yym(p,:)=uyyy';
end
pcolor(Yym,Xxm,rupall);
hold on;
plot([0 0],[0 50],'white')
plot([0 50],[0 50],'white')
xlabel('distance along the fault');
ylabel('seperation between two sources');
% cd ..
shading flat;
save (['2source_rupall'],'rupall');
saveas(gcf102,['2source_sameGR_x' ],'fig');
saveas(gcf102,['2source_sameGR_x' ],'eps');

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