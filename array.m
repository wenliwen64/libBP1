% clear all;
% close all;
function widthA = rour(P,sm,rm,h,w,N)
%sqr = @(x) x.^2;
%w=2*pi*[1 2 3 4 5];
%N=100;
%h=0.2;
P=[ 1 0 0 0 0];
s=linspace(0,0,sm,N);
r=linspace(0,rm,N);
A=zeros(1,length(s));
rou=zeros(1,length(r));
flag=1;
widthT=0;
for j=1:length(r)
     rou(j)=sum(P.*r(j).^((0:length(P)-1)));
            
end


  figure(12);
  hold all;
  %hold on;
for k=1:length(w)
    if flag==0
        break;
    end
    for i=1:length(s)
        
        for j=1:length(r)
            
            A(i)=A(i)+rou(j)*exp(1i*r(j)*s(i)*w(k));
        end
    end
    AA=abs(A).^2;
  
    
   plot(s,AA);
    
    width=interp1(AA,s,(max(AA)-min(AA))*h+min(AA));
    for i=1:length(s)
        if s(i)>width&&A(i)>(max(AA)-min(AA))*h+min(AA)
            flag=0;
        end
    end
   widthT=widthT+width;
end


if flag==0
    widthA=0;
else
    widthA=widthT/length(w);
end