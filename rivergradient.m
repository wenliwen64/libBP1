clear;
figure;
subplot(1,1,1);  
distance = [15 13 16 17 11 13 15 15 14 18 13 11 13 25 16 18 13 18 13 16 14 14 15 13 12];
distance1 = [13 12 13 15 2 1 6 14 7 8 6 7 12 11 12 10 6 7 4 4 5 6 3 13 13]
x(1)=0;
for i = 1:1:25,
    x(i+1)=distance(i)+x(i),
end
data=[800:-40:-200];
plot(x*1000/40,data,'RED')   

y(1)=0;
for i = 1:1:25,
    y(i+1)=distance1(i)+y(i),
end
[p,s] = polyfit(x,data,1)
[p1,s1] = polyfit(y,data,1)
 
for i = 1:1:25,
    d1(i)=(40*0.348)*40/distance1(i)/1000,
end

for i = 1:1:24,
    d2(i)=(d1(i+1)-d1(i))/distance1(i)/1000,
end
data
d2
y
p
p1
%plot(y(3:26)*1000/40,d2,'BLUE');


plot(x*1000/40,data,'RED',y*1000/40,data,'BLUE',x*1000/40,polyval(p,x),'BLACK',x*1000/40,polyval(p1,x),'BLACK')   
xlabel('Distance/meters');
ylabel('second derivitive');
title('Rurnace Creek second derivitive profile');

x

%f = fittype('a*x+b','problem','n','options',s);

%plot(x,data,'*',x,polyval(a,x),'o')