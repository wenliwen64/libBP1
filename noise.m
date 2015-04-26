clear;
load 'a'105;
load 'a'117;
aa1=a105(:,2);
aa2=a117(:,2);

meanaa1=sum(aa1)/length(aa1);
a1=aa1-meanaa1;
meanaa2=sum(aa2)/length(aa2);
a2=aa2-meanaa2;

z1=1.07*exp(i*2*pi*0.135/1);
cz1=conj(z1);
z0=1.08*exp(i*2*pi*0.135/1);
cz0=conj(z0);
B=[z0*cz0 -(z0+cz0) 1];
A=[z1*cz1 -(z1+cz1) 1];
a11=filter(B,A,a1)-a1;
a22=filter(B,A,a2)-a2;
[xco, lags] = xcorr(a11,a22);
plot(lags, xco);
figure;
subplot(2,1,1);
plot(a105(:,1),a1);
bbb=[lags' xco];

xlim([0 2000]);
subplot(2,1,2)
plot(a105(:,1),a11);
xlim([0 2000]);