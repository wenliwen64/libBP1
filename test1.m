% test the radiation pattern
close all

clear;
t = linspace(0,360,100);


     x=RadiationPattern(0,90,0,180-45,t);
figure(1)
       %polar(t/180*pi,(x(:,2).^2+x(:,3).^2).^0.5','*');
         polar(t/180*pi,abs(x(:,1))','*');
         
         
%          figure(2)
%          hold on;
% % U = RadiationPatte%        
% for i=1:100
%     if x(i,1)>0
%         polar(t(i)/180*pi,x(i,1)','*');
%     else
%         polar(t(i)/180*pi,x(i,1)','r*');
%     end
% end
