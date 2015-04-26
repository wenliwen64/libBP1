close all;
clear all;
station=['FUNV';'CRUV';'BAUV';'VIRV';'DABV';'MONV';'SIQV';'SANV';'QARV';'CURV';'CAPV';'SOCV';'CUPV';'JACV';'RIOV';'TURV';'TERV';'ELOV';'LUEV';'PCRV';'PRGV';'MAPV'];
 lat=[10.4692   10.6163    8.9433   10.5028   10.9218   11.9550   10.6488    9.5008   10.2065   10.0130    7.8647    8.2842   10.0563   11.0872    8.0618   10.4495    9.9637    7.0010 5.8432   10.1633    8.7600    9.8308];
 lon=-[66.8102   63.1842   68.0415   72.4060   70.6362   69.9703   69.8078   69.5363   70.5237   69.9612   72.3143   70.8563   65.7877   68.8298   61.8170   67.8395   69.1917   69.4833 61.4612   64.5897   64.6455   68.457];

%load /home/lsmeng/matlab/haitiUS/gf1/isap91_13/phtimefk;
load ptimes;
ttime=tt;
clear tt;
lat0=18.45;
lon0=-72.45;

sr=100;
load haiticoast.dat;
% main_normal.eps
%%%%%%%%%%%%%%%% for main event
load venx;
B=1:22;

    %%%%%%%%%%%%% M6 aftershock
% load retven60;
% xs=retven60.x;
% lat=retven60.lat;
% lon=retven60.lon;
% station=retven60.sta;
% B=[1:8 10:13 15 16];
%%%%%%%%%%%%



 B=[1 2 3 5 6 7 8 10 13:19 22];
nel=length(B);
x0=xs(B,:);
r=zeros(nel,2);
r(:,1)=lon(B);
r(:,2)=lat(B);
rdis=distance11(r(:,2),r(:,1),lat0,lon0,1)*180/pi;
nm=station(B,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BB,AA]=butter(4,[0.20 1]/(sr/2));
for i=1:nel
    
   x0(i,:)=filter(BB,AA,x0(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%


    for i=1:nel
          x01(i,:)=decimate(x0(i,:),20);
          
    end
    clear x0;
    x0=x01;
    sr=sr/20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nel
% x0(i,:)=x0(i,:)/std(x0(i,180*sr:185*sr));
% end
  
x01=x0;
%      load /home/lsmeng/matlab/chile/figfirday/y2ven.mat;
%     x01=zeros(nel,length(x0));
%     for i=1:nel
%           x01(i,:)=specshift(x0(i,:),y2(i)*sr);
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t3=(0:length(x0)-1)/sr;
for i=1:nel
    time(i,:)=t3;
end 
%%%%%%%%%%%

  Opr=readAllSac();
    Opr.lat0=18.45;
    Opr.lon0=-72.45;
    Opr.sr=sr;
  ret=struct('r',r,...%the position of the stations
    'nm',char(nm),...%name of the stations
    'rdis',rdis,...%epicentral distance
    'time',time,...% time matrix
    't1',0,...% theoretical arrival times
    'x',x0,...% original data 
    'x1',x0,...
    'sr',sr,...% sample rate
    'ori',180,...% origin time
    'opr',Opr);% attached the option object to read parameters.






%ret=seismoAlign(ret);




figure(21);

hold on;
for jj=1:nel
%plot(t3,x01(jj,:)/std(x01(jj,280*sr:320*sr))/8+(jj));
plot(t3,x0(jj,:)/std(x0(jj,180*sr:185*sr))/8+(jj),'r');
%plot(t3,x0(jj,:)/max(x0(jj,:))/8+(jj),'r');
end

xlim([170 250]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(22);
hold on;
for jj=1:nel
plot(t3,x01(jj,:)/std(x01(jj,180*sr:190*sr))/8+(jj));
%plot(t3,x0(jj,:)/std(x0(jj,280*sr:320*sr))/8+(jj),'r');
end
xlim([170 250]);
    clear x0;
    x0=x01;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [nel len ]=size(x0);
   time=zeros(nel,len);
    for i=1:nel
        time(i,:)=t3;
    end
    
    


    
    
    
    
    
    
    
    
    Opr=readAllSac();
    Opr.lat0=18.45;
    Opr.lon0=-72.45;
    Opr.sr=sr;
  ret=struct('r',r,...%the position of the stations
    'nm',char(nm),...%name of the stations
    'rdis',rdis,...%epicentral distance
    'time',time,...% time matrix
    't1',0,...% theoretical arrival times
    'x',x0,...% original data 
    'x1',x0,...
    'sr',sr,...% sample rate
    'ori',180,...% origin time
    'opr',Opr);% attached the option object to read parameters.

plotopr=plotAll();
plotopr.normal=[10 15 8];
plotopr.view=[0 50];
plotopr.plotAlignbool=false;

%plotAll(ret,plotopr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotSpec(ret,10,50,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=0;
doaopr=doaAll1();
doaopr.method='music';
doaopr.freqBand=[0.4 0.4 0.4]*3/4;
doaopr.ps=40;
doaopr.qs=40;
doaopr.step=1;
doaopr.begin=-60;
doaopr.over=-15;
doaopr.windowLength=60;
doaopr.slidingWindowViewRange=[20 20 1];
doaopr.uxRange=[-1.5 0.9];
doaopr.uyRange=[-1.2 1.2];
doaopr.projectionRange=[-20 100];
doaopr.constantBias=0;
%doaopr.beamOrientation=-(18.5558-18.4888)/(72.4744-72.4313);
doaopr.beamOrientation=125;
doaopr.saveProj=NN;
doaopr.Mm=2;
doaopr.saveDir='/home/lsmeng/matlab/haitiUS/gf/save_half/';

doaAll1(ret,doaopr);


















