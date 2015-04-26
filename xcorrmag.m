clear all;
close all;
load m5data;
m1=m5-5;
clear m5;
load m5the;
m2= m5;
figure(2);
plot(1:length(m5),m1,1:length(m5),m2)
C=xcorr(m1,m2,'coeff');
figure(1);
plot(-length(m5)+1:length(m5)-1,C);
