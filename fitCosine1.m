clear all;
close all;
d=textread('fletcher_etal_2006_fig6_events.txt', '' , 'headerlines', 1); 
figure(1);
baz=d(:,8);
tbaz=d(:,19);
residue=baz-tbaz;
plot(tbaz,residue,'r.');
xlim([0 360]);
az=0:360;

A=1:60;
T=1:360;
P=zeros(length(A),length(T));
for i=1:length(A)
    for j=1:length(T)
%         P(i,j)=sum((residue-A(i)*cosd(tbaz+T(j))).^2);
        P(i,j)=sum((residue-A(i)*cosd(baz+T(j))).^2);

    end
end

ind=peakfit2d(-P);
A0=interp1(1:length(A),A,ind(1));
T0=interp1(1:length(T),T,ind(2));

hold on;
plot(az,A0*cosd(az+T0));
% xlabel('True back azimuth (^o )');
xlabel('observed back azimuth (^o )');

ylabel('Residue (observed - true back azimuth) (^o )');
legend('Residue','Best-fitting consine curve');
title(['A0 = ' num2str(A0) ' T0 = ' num2str(T0)]);
% figure(2);
% pcolor(-P);