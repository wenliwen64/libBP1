function [uxx uyy val distpeak max1]=tprfslow(ux,uy,Pm,bux,buy,lat0,lon0,latm,lonm,dist,faultorientation) 
%   dist0=distance22(lat0,lon0,latm,lonm)*110.49;
% a0=faultorientation-azimuth(lat0,lon0,latm,lonm);
% dist00=dist0*tand(a0);
n=length(dist);
az=zeros(1,n);
val=zeros(1,n);
slowpeak=min([ sqrt(bux^2+buy^2) 0.3]);
%slowpeak=sqrt(bux^2+buy^2)

% slowpeak=0.4;

dlat=dist*cosd(faultorientation)/110.49;
dlon=dist*sind(faultorientation)/cosd(latm)/110.49;
for i=1:n
%     az(i)=faultorientation-180/pi*atan2((dist0*sind(a0)),(dist0*cosd(a0)-dist(i)));
%      az(i)=faultorientation-90+atan()
%     dista=sqrt(dist0^2+dist(i)^2-2*dist0*dist(i));
%     a00=real(asind(sind(a0)*dist0/dista));
%     az(i)=faultorientation-180+a00;
    
    az(i)=azimuth(lat0+dlat(i),lon0+dlon(i),latm,lonm);
    uxx(i)=slowpeak*sind(az(i));
    uyy(i)=slowpeak*cosd(az(i));
    val(i)=interp2(ux,uy,Pm',uxx(i),uyy(i));
   
end
max1=max(val);
distpeak=dist(val==max1);
% figure(20);
% plot(lon0,lat0,'k*',lonm,latm,'r.',lon0+dlon(1),lat0+dlat(1),'g*',lon0+dlon(end),lat0+dlat(end),'r*');
% 111