figure(1);
hold on;
a=[1:11]*8.3;
b=[8 17 26 33 42 51 61 68 79 89 96]
% b=[7 11 18 23 28 31 35 39 43 47 52 55];
% d=[7 10 19 24 26 30 35 38 41 45 53 55];
plot(a,b,'*-');
%plot(a,d,'g*-');
c=polyfit(a,b,1);

plot(a,a*c(1)+c(2),'r.-');