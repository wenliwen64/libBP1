load aria01wt11;
A=aria01wt11;
lat=A(:,1);
re=A(:,4);
re1=A(:,5);
figure(3);
plot(lat,re-mean(re),lat,re1-mean(re1),'r');
