close all;
close all;
for NN=1:12
    close all;
cd '/home/lsmeng/matlab/haitiUS/gf'
system(['makesynm.sh ' num2str(NN) ' ' num2str(NN*3)]);
opr=readAllSac();
opr.FileList='filelist';
opr.lat0=18.45;
opr.lon0=-72.45;
opr.sr=5;
opr.ori=100;
%opr.align.timeShiftKnown=load ('/home/lsmeng/matlab/haitiUS/gf/syn/arrival000')-200;
opr.align.shiftT1bool=true;
opr.alignbool=true;
opr.bpbool=true;
opr.bp=[0.5 1 ];
ret=readAllSac('/home/lsmeng/matlab/haitiUS/gf/syn/',opr);
plotSta(ret);

plotopr=plotAll();
plotopr.normal=[20 20 2];
plotopr.view=[-100 500];
plotopr.plotAlignbool=true;

plotAll(ret,plotopr);

doaopr=doaAll();
doaopr.method='stacking';
doaopr.freqBand=[1 1];
doaopr.ps=30;
doaopr.qs=45;
doaopr.step=1;
doaopr.begin=101;
doaopr.over=160;
doaopr.windowLength=5;
doaopr.slidingWindowViewRange=[20 20 0.25];
doaopr.uxRange=[-1.2 1.2];
doaopr.uyRange=[-1.8 1.8];
doaopr.projectionRange=[-20 100];
%doaopr.beamOrientation=-(18.5558-18.4888)/(72.4744-72.4313);
doaopr.beamOrientation=-(18.4086-18.243)/(72.4084-72.3259);
doaopr.saveProj=NN;
doaopr.saveDir='/home/lsmeng/matlab/haitiUS/gf/save2/';
doaAll(ret,doaopr);
end
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