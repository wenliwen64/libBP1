clear all;
close all;
% cd /home/lsmeng/matlab/bpsynthetic/;
%addpath /u/home/w/wenliwen/libBP1
%addpath /Users/lsmeng/Dropbox/matlab/formalfunc%this one needs to be up than libBP
opr=readAllSac();
opr.bpbool=false;
opr.lon0=-70.817;
opr.lat0=-19.642;
opr.dep0=20;
opr.bp=[0.01 1];
opr.snrFilterbool=false;
opr.ori=300;
opr.snrFilter=[0.1 0.5 2 -20 -10 100 130];
opr.sr=10;
load /u/scratch/w/wenliwen/Data/2014_04_01_1_139/2014_04_01_1_139.mat 
ret=new_struct;
ret.opr=opr;
ret.timeshiftall=zeros(1,length(ret.lon));
ret.lon=round(ret.lon*1e4)/1e4;
% ret.xori=ret.x;
% ret=rmfield(ret,'x');
tmp=load(' /u/home/w/wenliwen/newFolder/libBP1/nchileCali.mat');
ret1=tmp.ret;
ret1.lon=round(ret1.lon*1e4)/1e4;
ret=matchshiftAll(ret,ret1);
% ret.x=ret.xori;%%
% load nchile1;
if isfield(ret,'time')
    ret=rmfield(ret,'time');
end
[n m]=size(ret.xori);
ret.lon0=ret.opr.lon0;
ret.lat0=ret.opr.lat0;
ret.dep0=ret.opr.dep0;
ret.parr= 20;
ret.begin= 0;
ret.end=round(m/ret.sr)-20;
ret.step= 1;
ret.ps= 40;
ret.qs= 40;
ret.lonrange= [-0.5 1];
ret.latrange= [-1 0.5];
ret.fl= 0.5;
ret.fh= 2;
ret.fs=2;
ret.win= 5;
ret.dirname= 'nchile_TA';
ret.Nw= 3;
fl=0.5;
fh=2;
save('nchileinput_1','ret','-v7.3');
% ret.x=ret.xori;
% [BB,AA]=butter(4,[fl fh]/(ret.sr/2));
% for i=1:n
%     
%      ret.x(i,:)=filter(BB,AA,ret.x(i,:));
% end
% 
% ret.scale=1.50;
% ret.lt=300;
% ret.ht=350;
% figure(12);
% plotAll1(ret);
% xlim([200 400]);
% runteleBPcont(ret);
