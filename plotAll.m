function ret=plotAll(ret,plotOpr)
% function ret=plotAll(ret,plotOpr)
% plot the data in two windows : original and aligned
def=struct('normal',[-20 20 8],...% normalization range and normalization ratio with respect by std
            'view',[-10 50],...% view range of the seismograms
            'bp',[0 0],...%filtering
            'plotAlignbool',false); % whether to plot the aligned data or not
if ~nargin %user wants default options, give it and stop
    ret = def;
    return
elseif nargin<2
    plotOpr=def;
end
x=ret.x;
x1=ret.x1;
[nel lengthSam]=size(ret.x);
sr=ret.sr;
ori=ret.ori;
 t3=(0:length(x)-1)/sr;
 
 
 
h=figure(20);
close(h);
figure(20);
hold on;

if plotOpr.bp(1)~=plotOpr.bp(2)
    [BB,AA]=butter(4,plotOpr.bp/(sr/2));
    for i=1:nel
        
        x(i,:)=filter(BB,AA,x(i,:));
        x1(i,:)=filter(BB,AA,x1(i,:));
    end
end

for jj=1:nel
    %plot(t3,x01(jj,:)/std(x01(jj,280*sr:320*sr))/8+(jj));
    if plotOpr.normal(1)==plotOpr.normal(2)
        plot(t3,x(jj,:)/std(x(jj,:))/plotOpr.normal(3)+(jj),'r');
    else
        plot(t3,x(jj,:)/std(x(jj,(ori+plotOpr.normal(1))*sr:(ori+plotOpr.normal(2))*sr))/plotOpr.normal(3)+(jj),'r');
    end
end
xlim(plotOpr.view+ori);
    
if plotOpr.plotAlignbool==true% plot the aligned data
    
        h=figure(21);
        %close(h);
        figure(21);
hold on;
for jj=1:nel
    jj
    size(t3)
    size(x1)
    if plotOpr.normal(1)==plotOpr.normal(2)
        plot(t3,x1(jj,:)/std(x1(jj,:))/plotOpr.normal(3)+(jj),'b');
    else
        plot(t3,x1(jj,:)/std(x1(jj,(ori+plotOpr.normal(1))*sr:(ori+plotOpr.normal(2))*sr))/plotOpr.normal(3)+(jj),'b');
    end
    %plot(t3,x0(jj,:)/std(x0(jj,280*sr:320*sr))/8+(jj),'r');
end
xlim(plotOpr.view+ori);
   
end