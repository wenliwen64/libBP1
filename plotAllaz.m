function plotAll1(ret)
x=ret.x;

sr=ret.sr;
[n m]=size(x);
t=(1:m)/sr;

hold on;
for i=1:n
    
     plot(t,x(i,1:m)/(std(x(1,ret.lt*sr:ret.ht*sr)))*ret.scale+ret.az(i));
%           plot(t,x(i,1:m)/(std(x(i,:)))*ret.scale+i*10);
%     text(t(1*sr),i*10,num2str(ret.rdis(i)),'FontSize',5,'color','black');
%      text(t(55*sr),i*10,num2str(ret.az(i)),'FontSize',5,'color','black');
%       text(t(52*sr),i*10,num2str(ret.rdis(i)),'FontSize',5,'color','black');
%      text(t(50*sr),i*10,num2str(ret.nm(i,:)),'FontSize',5,'color','black');
%     text(t(70*sr),i*10,num2str(1e4*std(x(i,1:m))),'FontSize',5,'color','black');
end