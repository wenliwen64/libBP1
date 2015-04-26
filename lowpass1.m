
close all;
clear all;
figure(1);
sr = 200; %sampling rate /sec

%t = linspace(0,1,sr);
%x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
dur=10;
nel=12;
x1=zeros(nel,dur*sr);
k=1;
filename='2760704k4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6);
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x1(k,:) pr]=load_bbdata([filename '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x1(k,:) pr]=load_bbdata([filename num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x1(k,:)=x1(k,:)-mean(x1(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end
x2=zeros(nel,dur*sr);
k=1;
filename1='2740502k4.p';
[xtext pr]=load_bbdata([filename1 '13' ],[],dur);

start=pr.t0;
start(6)=start(6);

for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x2(k,:) pr]=load_bbdata([filename1 '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x2(k,:) pr]=load_bbdata([filename1 num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x2(k,:)=x2(k,:)-mean(x2(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end
x=zeros(nel,dur*sr);
x=x2;
%x=xx(1001:1500);
t=linspace(0,dur,dur*sr);
%subplot(3,2,1);                   
%[B,A]= cheby1(3,0.9, [16.9 23.7]/sr); %construct the filter for 8 to 12 hz
figure(26);

D = fdesign.lowpass('Fp,Fst,Ap,Ast',0.2, 0.22, 1, 60);
    butterbps = design(D,'equiripple');
%bps = fdesign.highpass(9.5, 10.8, 11.2,  .7, 7, 200);
%butterbps=butter(bps);
figure(27);
[B,A]= butter(2, [12.5 14.5]/sr); %construct the filter for 8 to 12 hz
[H,F]=freqz(butterbps,1024,sr);
plot(F,abs(H));
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('Chebyshev filter 4th order Frequency response');

%subplot(3,2,2);                  
%[H1,F1]=impz(B,A,1024,sr);
%plot(F1,abs(H1));
%xlabel('Time(s)');
%ylabel('Amplitude');
%title('Chebyshev filter 4th order impulse filter');

%subplot(3,2,3);   
y=zeros(nel,dur*sr);
for i=1:nel
    y(i,:)= filter(butterbps,x(i,:));
end
y1 = filter(butterbps,x(1,:));
y2= filter(B,A,x(2,:));
y3= filter(B,A,x(3,:));
%plot(t,x,t,y,'r');
%xlabel('Time(s)');
%ylabel('Amplitude');
%title('signals time series');
figure(2);
%subplot(3,2,4);  
a=fft(x(1,:));
b=fft(y1);
tt=linspace(0,sr,sr*dur);
plot(tt,abs(a),tt,abs(b),'r','markersize',1);
xlim([0 0.5*sr]);
xlabel('Frequency(Hz)');
ylabel('Amplitude');
title('signals in frequency domain');
figure(3);
for i=1:nel
subplot(nel,1,i)
plot(t,y(i,:),'b');
end
% subplot(3,1,2);
% plot(t,x(2,:),'b');
% subplot(3,1,3);1
% plot(t,x(2,:),'b');
%xlim([1 2]);
%subplot(3,2,5);                  
%zeros=roots(B);
%poles=roots(A);
%zplane(zeros,poles);
%xlabel('Real part');
%ylabel('imaginary');
%title('zeros and poles');
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;


sd274=[ -0.0214   -0.0219    0.0371   -0.0034   -0.0012   -0.0010  0.0133 -0.0110   -0.0220    0.0020    0.0333   -0.0040];%274
  sd291=[-0.0152   -0.0227    0.0357   -0.0092   -0.0072    0.0002    0.0157  -0.0074   -0.0212    0.0028    0.0338   -0.0053] ;%291

 sd282=[   -0.0187   -0.0213    0.0402   -0.0072   -0.0085   -0.0010 0.0145   -0.0075   -0.0179   -0.0008    0.0304   -0.0022];% 282
  %sd=[     0.0035    0.0000    0.0039   -0.0079   -0.0032   -0.0074    0.0057    0.0039    0.0010   -0.0003    0.0013   -0.0006];% 282
 %sd=(sd274+sd291+sd282)/3;
     %sd= [0.0915    0.0648    0.1022    0.0066         0    0.0048   -0.0332   -0.0740   -0.1245   -0.0673   -0.1216   -0.1729]
%sd=[    0.1017    0.0776    0.0535    0.0054         0   -0.0027   -0.0562   -0.0749   -0.1151   -0.0749   -0.1605   -0.1793]%y1 282
%sd=[    0.0884    0.0576    0.0954    0.0033         0    0.0053   -0.0243   -0.0747   -0.1395   -0.0751   -0.1331   -0.1630]%282 uz=0
  sd=zeros(1,12);
  ux=0.5;
  uy=-0.6;
  %sd=[ux uy (1/0.7^2-ux^2-uy^2)^0.5]*r';
  
  ux1=0.5;
  uy1=-0.3;
  %sd1=[ux1 uy1 (1/0.7^2-ux1^2-uy1^2)^0.5]*r';
  
window=0.5;
step=0.05;
k=1;
ts1=2.0;
figure(25);
for i=1:nel
subplot(nel,2,2*i-1);
plot(t(round((ts1)*sr)+1:round((ts1)*sr)+window*sr),x(i,round((ts1+sd(i))*sr)+1:round((ts1+sd(i))*sr)+window*sr));

end

for i=1:nel
subplot(nel,2,2*i);
plot(t(round((ts1)*sr)+1:round((ts1)*sr)+window*sr),y(i,round((ts1+sd(i))*sr)+1:round((ts1+sd(i))*sr)+window*sr));

end

for ts1=1.5:step:2.6
%ts1=1;
%for ts1=0:step:dur-window
clear yy;
ts2=ts1+window;
yy=zeros(nel,window*sr);
for j=1:nel
    
yy(j,1:window*sr)=x(j,round((ts1+sd(j))*sr)+1:round((ts1+sd(j))*sr)+window*sr);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
[X,Uv,A]=musicDOAdft(yy,window,ts1);
end
%[Pm,U,A]=musicDOAslowness(yy,window);

%end
% theta1(k)=the(1);
% phi1(k)=ph(1);
% std_all=0;
% for i=1:nel
% std_all=std_all+std(yy(i,:))^2;
% end
% Em(k)=Emm/std_all/window/sr;
% k=k+1;
% end
% 
% 
% tss=linspace(0,dur-window,k-1);
% figure(4);
% N=5;
% subplot(N,1,1);
% plot(tss,theta1*180/pi);
% xlabel('time(sec)');
% ylabel('theta(degree)');
% subplot(N,1,2);
% plot(tss,phi1*180/pi);
% xlabel('time(sec)');
% ylabel('phi(degree)');
% subplot(N,1,3);
% plot(tss,Em);
% subplot(N,1,4);
% plot(t,x(5,:));
% subplot(N,1,5);
% plot(t,y(5,:));


