function [dist val x yy]=tprf(Ym,Xm,Pm,lonm,latm,lat0,lon0,a,b,pji)
%function : [dist val x yy]=tprf(Ym,Xm,Pm,lonm,latm,lat0,lon0,a,b,pji)
%to project the peak of matrix on the faultline
%input:
% Ym,Xm,Pm :the backprojection matrix.
% lonm, latm :the position of the peak in Pm
% lat0 lon0 :epicenter
% a beam :orientation
% b fault :orientation
% pji the :range along the fault relative to the epicenter
%output:
% dist: distance vecter along the fault
% val: interpolated values along the fault
% x yy: the position vectors(lon and lat ) of the fault plane 
fl=(pji(2)-pji(1))/111.94;
latl=fl*sind(b);
lonl=fl*cosd(b)/cosd(lat0);
% latl=0.0887;
% lonl=0.2512;
% lat0=18.44;
% lon0=-72.57;
% latm=18.55;
% lonm=-72.52;
lati=linspace(latm-latl*4,latm+latl*4,800);
loni=linspace(lonm-lonl*4,lonm+lonl*4,800);
val0=interp2(Ym,Xm,real(Pm),loni,lati,'linear',0);
for i=1:length(lati)
    x1=lon0;
    y1=lat0;
%     x2=lon0-0.2512*2;
%     y2=lat0-0.0887*2;
    x2=lon0-lonl*2;
    y2=lat0-latl*2;
    x3=lon0-lonl*10;
    y3=lat0-latl*10;
    %a=-(18.5558-18.4888)/(72.4744-72.4313);
    x0(i)=loni(i);
    y0(i)=lati(i);
    a1=(y2-y1)/(x2-x1);
    b1=(x1*y2-x2*y1)/(x1-x2);
    a2=a;
    b2=y0(i)-a*x0(i);
    x(i)=(b2-b1)/(a1-a2);
    yy(i)=(a2*b1-a1*b2)/(a2-a1);
    
    dist0(i)=(-distance22(yy(i),x(i),y3,x3)+distance22(y1,x1,y3,x3))*111.1944;
end
dist=pji(1):pji(2);
val=interp1(dist0,val0,dist,'linear',0);

%plot(x1,y1,'r*',x3,y3,'r*',x0,y0,'b',x,yy,'black',lonm,latm,'g*');

