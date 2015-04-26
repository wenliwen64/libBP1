function plotAll1az(ret)
x=ret.x;

sr=ret.sr;
[n m]=size(x);
t=(1:m)/sr;

hold on;
if isfield(ret,'ploty')
    ploty=ret.ploty;
else
    ploty=1:n;
end
if isfield(ret,'plotcolor');
    plotcolor=ret.plotcolor;
end
for i=1:n
    if isfield(ret,'plotcolor');
        if ret.plotcolor(i)==0
            plotcolor='yellow';
            
        elseif ret.plotcolor(i)==-1
            plotcolor='red';
            
        elseif ret.plotcolor(i)==1
            plotcolor='blue';
            
        end
    else
        
        plotcolor='b';
    end
    plot(t,x(i,1:m)/(std(x(1,ret.lt*sr:ret.ht*sr)))*ret.scale+ploty(i),plotcolor);
    %           plot(t,x(i,1:m)/(std(x(i,:)))*ret.scale+i*10);
    %     text(t(1*sr),i*10,num2str(ret.rdis(i)),'FontSize',5,'color','black');
    %      text(t(55*sr),i*10,num2str(ret.az(i)),'FontSize',5,'color','black');
    %       text(t(52*sr),i*10,num2str(ret.rdis(i)),'FontSize',5,'color','black');
    %  text(t(450*sr),ploty(i),num2str(ret.rdis(i)),'FontSize',5,'color','black');
    %     text(t(70*sr),i*10,num2str(1e4*std(x(i,1:m))),'FontSize',5,'color','black');
end