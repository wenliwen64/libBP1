function [dist val x yy]=tpr(Ym,Xm,Pm,lonm,latm,a)

lat0=18.44;
lon0=-72.57;
% latm=18.55;
% lonm=-72.52;
lati=linspace(latm-0.0887*4,latm+0.0887*4,400);
loni=linspace(lonm-0.2512*4,lonm+0.2512*4,400);
val0=interp2(Ym,Xm,real(Pm),loni,lati,'linear',0);
for i=1:length(lati)
    x1=lon0;
    y1=lat0;
    x2=lon0-0.2512*2;
    y2=lat0-0.0887*2;
    x3=lon0-0.2512*10;
    y3=lat0-0.0887*10;
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
dist=0:100;
val=interp1(dist0,val0,dist,'linear',0);
plot(x1,y1,'r*',x2,y2,'g*',x3,y3,'r*',x0,y0,'b',x,yy,'black',lonm,latm,'g*')