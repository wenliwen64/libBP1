function [dist val lon lat]=tprfa(Ym,Xm,Pm,latm,lonm,lat0,lon0,a,b,pji,N)
% b



dist=linspace(pji(1),pji(2),N);
val=zeros(1,length(dist));
lon=zeros(1,length(dist));
lat=zeros(1,length(dist));

 latj=4*sind(a);
 lonj=4*cosd(a)/cosd(lat0);
for i=1:N
    
    lat(i)=lat0-dist(i)*sind(b)/111.94;
    lon(i)=lon0-dist(i)*cosd(b)/cosd(lat0)/111.94;
   
    latk=linspace(lat(i)-latj,lat(i)+latj,100);
%     lonk=linspace(lon(i)-lonj,lon(i)+lonj,100);
%     val0=interp2(Ym,Xm,real(Pm),lonk,latk,'linear',0);
%    
%     valmax=max(val0);
%     for j=1:length(latk)
%         if val0(j)==valmax
%             lat1=latk(j);
%             lon1=lonk(j);
%         end
%     end
%     
%     lat2=linspace(lat1-0.5*sind(b),lat1+1*sind(b),100);
%     lon2=linspace(lon1-0.5*cosd(b)/cosd(lat0),lon1+1*cosd(b)/cosd(lat0));
%     val1=interp2(Ym,Xm,real(Pm),lon2,lat2,'linear',0);
%     %val(i)=max(val0)^2/max(val1);
%     val(i)=max(val0);
end
     val=interp2(Ym,Xm,real(Pm),lon,lat,'linear',0);

  latm
  lonm
  latl=linspace(latm-latj,latm+latj,100);
  lonl=linspace(lonm-lonj,lonm+lonj,100);
  plot(lonl,latl,'white','LineWidth',1);