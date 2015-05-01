clear all;
close all;
% cd /home/lsmeng/matlab/bpsynthetic/;
opr=readAllSac();
opr.bpbool=false;
opr.lon0=-70.817;
opr.lat0=-19.642;
opr.bp=[0.01 1];
opr.snrFilterbool=false;
opr.ori=300;
opr.snrFilter=[0.1 0.5 2 -20 -10 100 130];
opr.sr=10;
% ret=readAllSac('/Users/lsmeng/wk/matlab/bp/nchile/data/',opr);

% 
%  save Mexico032012;
% load Mexico032012_1;
% load el5;
load 201504251500_240.v1.mat;
%   load can2;
%  ret.xori=ret.x;
 ret.x=ret.xori;
% load kyushu.mat;
% load tohoku.mat;
%  ret.xori=ret.x;
% save('nchile','ret');
%  load HonshuORF1;
 load ptimes;
% I=find(ret.lon<139.5);
% ret=orderAll(ret,ret.rdis,I);
% I=find(ret.lon>136.6);
% ret=orderAll(ret,ret.rdis,I);
%  I=find(ret.lat<36.5);
%  ret=orderAll(ret,ret.rdis,I);
%   I=find(ret.lat>34.5);
%  ret=orderAll(ret,ret.rdis,I);
% I=find(ret.az>180.0);
% ret=orderAll(ret,ret.rdis,I);

 [n m]=size(ret.x);
% % I=find(ret.az>312);
% % ret=orderAll(ret,ret.rdis,I);
% % I=find(ret.az<320);
% % ret=orderAll(ret,ret.az,I);
% I=50+([1:62 64 66:74 76:81 83:166]);
% I=51:216;
% I=[1:3 5:61 63:76 78:n];


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ret=orderAll(ret,ret.rdis);
  I=2:8;
%  ret=orderAll(ret,0,I);
  [n m]=size(ret.x);
  figure(1);
plotSta(ret);
%  
%    x1=zeros(n,length(200*ret.sr:1600*ret.sr));
%    for i=1:n
%        x1(i,:)=ret.x(i,200*ret.sr:1600*ret.sr);
%    end
%        ret.x=x1;
% %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %        [n m]=size(ret.x);
% %      
% %        []
% %    ts=(1:m)/ret.sr;
% %  pcolor(ts,1:n,ret.x);   
% %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ptimesSocalP
% shiftP=interp1(rr,tt,ret.rdis);
% load ptimesSocalS
% shiftS=interp1(rr,tt,ret.rdis);
% shift=shiftS-shiftP;
            for i=1:n
%               ret.x(i,:)=specshift(ret.x(i,:),ret.sr*(shift(i)-shift(1)));
%                          ret.x(i,:)=specshift(ret.x(i,:),ret.sr*(-100));
%            ret.x(i,:)=specshift(ret.x(i,:),ret.sr*((ret.t1(i)-ret.t1(1))*3600*24));
%            ret.x(i,:)=specshift(ret.x(i,:),-ret.sr*((ret.t1(i)-ret.t1(1))));

            end
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d
        ret.xori=ret.x;
% %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        load HonshuEU6;
% ret=shiftAll(ret,-40)
%      ret.x=ret.xori;
%      I=1:60;
%      ret=orderAll(ret,0,I);
%  load HonshuORF6;
% 
% %  I=find(ret.mxcr1>0.8);
% % ret=orderAll(ret,ret.rdis,I);
%  I=find(ret.az<47.3);
% ret=orderAll(ret,ret.rdis,I);
% 
% I=find(ret.az>28.0);
% ret=orderAll(ret,ret.rdis,I);
%  ret.x=ret.xori;
%   figure(15);
%  plotSpec(ret,400,700,1);
%  ret=orderAll(ret,ret.az);
%     [n m]=size(ret.x);
%     
%     %%%%%%%%%%%%%%%%%%%%%%

fl=0.25;
fh=1;
[BB,AA]=butter(4,[fl fh]/(ret.sr/2));
for i=1:n
    
     ret.x(i,:)=filter(BB,AA,ret.x(i,:));
end
for i=1:n
  
     ret.x(i,:)=ret.x(i,:)/std(ret.x(i,10*ret.sr:400*ret.sr));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a=[];
% for i=1:n
%       dd=std(ret.x(i,450*ret.sr:480*ret.sr))/std(ret.x(i,500*ret.sr:520*ret.sr));
%       cc=std(ret.x(i,500*ret.sr:520*ret.sr))/std(ret.x(i,450*ret.sr:480*ret.sr));
%     if dd< 0.3 && cc >0.5
%         a=[a i];
%     end
%         
% end
% 
% ret=orderAll(ret,0,a);
% [n m]=size(ret.x);
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% %       I=[1:200];
% % ret=orderAll(ret,0,I);
% % [n m]=size(ret.x);
% %%%%%%%%%%%%%%%%%%%%%%%%
% % 
align=seismoAlign();
align.ts11=300;
align.win=10;
align.lt=400;
align.range=5;
align.refst=145;
%    align.cutoff=0.8;
ret=seismoAlign(ret,align);
% save 2012Mexico_1
align=seismoAlign();
align.ts11=398;
align.win=10;
align.lt=100;
align.range=1;
align.refst=0;
   align.cutoff=0.8;
%   ret=seismoAlign(ret,align);
% ret=seismoAlign(ret,align);
% 
%   ret=seismoAlignMulit2(ret,align);
% % 
% % 
% % 
% % ret=seismoAlignMulit(ret,align);
% %  save SIEDCARafter1;
% 
% 
% ret=pickAll(ret,0.5,15);
% ret=pickAll(ret,0.5,30);

ret.scale=1.50;
ret.lt=300;
ret.ht=350;
% h=figure(10);
% % title(['aftershock' num2str(fl) '-' num2str(fh) 'Hz'] )
% % subplot(3,2,1);
figure(10);
plotAll1(ret);
xlim([-100 100]);
% ret=pickAll(ret,1,13);
% 
% % saveas(h,['aftershock_data' num2str(fl) '_' num2str(fh) 'Hz.eps'],'epsc');
% % saveas(h,['aftershock_data' num2str(fl) '_' num2str(fh) 'Hz.fig'],'fig');
% % h=figure(11);
% % for i=1:6
% %     t0=60+i*10;
% %     subplot(3,2,i);
% % corrmap(ret.x(:,(t0)*ret.sr:(t0+20)*ret.sr)',ret.nm(:,3:4));
% % if i~=1
% %  colorbar off;
% % end
% % shading flat;
% % xlabel([num2str(t0) '-' num2str(t0+20) 's' ]);
% % ylabel('');
% % title(['aftershock' num2str(fl) '-' num2str(fh) 'Hz']);
% % caxis([0.5 1]);
% % end
% % saveas(h,['aftershock_corr' num2str(fl) '_' num2str(fh) 'Hz.eps'],'epsc');
% % saveas(h,['aftershock_corr' num2str(fl) '_' num2str(fh) 'Hz.fig'],'fig');
% %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  save USleft2;
% % save HonshuEU5
 save('nepal_align1','ret');