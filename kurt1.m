function [skew val]=kurt1(x0,k0,type,flag,C)
if nargin == 1
    k0=0;
    type=1;
    flag=0;
elseif nargin == 2
    type=1;
    flag=0;
elseif nargin == 3
    flag=0;
end
x0=x0-mean(x0);
n=181;
theta=linspace(-90,90,n)+k0;
kurtosis=zeros(1,n);
m=length(x0);
% type=3;
for i=1:n
    
    x1=hilshift(x0,theta(i));
%     x1=detrend(x1);
%      xx1=abs(x1);
xx1=x1;
%     xx1=[x1' fliplr(x1')];
%     size(xx1)
%         xx1=xx1(round(m/20):end-round(m/20));
%        xx2=fliplr(xx1);
%     xx1=xx1;
    if type==1
%              kurtosis(i)=std(abs((diff(xx1)))) ;
        kurtosis(i)=mean(xx1.^4)/(mean(xx1.^2))^2-3; %kurtosis
    elseif type==2
%              kurtosis(i)=mean(abs(diff(xx1))) ;
          Ip=find(xx1>=0);
          In=find(xx1<0);
        kurtosis(i)=std(xx1(Ip))/mean(xx1(Ip))+std(xx1(In))/mean(-xx1(In));% std(abs)/mean(abs)
    elseif type==3
        kurtosis(i)=std((diff(xx1))) ;% mean(abs)
    end
    
end
% figure
[val index]=peakfit1d(-kurtosis);
val=-val
skew=interp1(1:length(theta),theta,index);
if flag==1
plot(theta,kurtosis,C);
% plot(skew,val,[C '.'])
end
disp(['skew0=' num2str(skew)]);
% if nargout==2
%     varargout{1}=val;
%     varargout{2}=skew;
% %     [val skew];
% elseif nargout ==1
%     varargout=skew;
% end
