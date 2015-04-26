clear;
figure;
title('plot of first-motion');
polar(0,0,'.')
hold on;
load mech_dat_1;
fm=mech_dat_1(:,1);
az=mech_dat_1(:,3);
ih=mech_dat_1(:,4);
for k=1:length(fm)
    
    if ih(k)>90
        az(k)=az(k)+180;
        ih(k)=180-ih(k);
    end
    az(k)=az(k)/180*pi;
    ri(k)=2^0.5*sind(ih(k)/2);
    if fm(k)==1
        polar(-az(k),ri(k),'*r');
    else
        polar(-az(k),ri(k),'ob');
    end
end