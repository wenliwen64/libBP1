figure(1);
%hold on;
I=[-24 -22 -21 -19 -17 -15 -13 -11 -9 -7 -5 -3 -1 0.5 2.5 4.5 7 9 11 13 15 16.5 18.5 20.5 22 24 25.5 27.5 29 30.5 32 33.5 35 36.5 38 39 40 42]*pi/180;
Ia=I(1:23);
Ib=I(1:23);
Ib=I(15:37);
lata=-15:7;
latb=0:22;
alpha=170*pi/180;
D=0;
for i=1:23
    ea(i)=180/pi*atan(tan(Ia(i))/sin(alpha-D));
    eb(i)=180/pi*atan(tan(Ib(i))/sin(alpha-D));
    theta(i)=180-ea(i)-eb(i);
end
for i=1:22
    
    ddl(i)=(theta(i+1)-theta(i))
end
plot(lata,theta,'r-*');
figure(2);
plot(lata(1:22),ddl,'-*');