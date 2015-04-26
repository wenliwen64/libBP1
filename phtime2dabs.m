function t=phtime2dabs(lat1,lon1,dep1,lats,lons,rr,dd,tt)
% f=1;
% lat1=lat0;
% lon1=lon0;
% lat2=ux(p);
% lon2=uy(q);
% tt=ttime;
% lats=lat;
% lons=lon;
r1=distance22(lat1,lon1,lats,lons);
t=interp2(rr,dd,tt,r1,dep1);


%ph=exp(-1i*2*pi*f*t);
%if max(r1)>90||max(r2)
% r=0:0.1:30;
% t=zeros(1,length(r));
% for i=1:length(r)
% rr=num2str(r(i));
% [status,result]=system(['phtimes 10 ' rr ' P']);
% t(i)=str2double(result(5:20));
% endfigu