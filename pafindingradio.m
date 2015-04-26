clear all;
close all;
sr = 200; %sampling rate /sec
%t = linspace(0,1,sr);
%x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
dur=12;
nel=12;
x=zeros(nel,dur*sr);
xe=zeros(nel,dur*sr);
xn=zeros(nel,dur*sr);
k=1;
t=linspace(0,dur,dur*sr);
fall=  ['2721715j';'2740502k';'3091728j';'3091358m';'2960421s';'2960752d';'2751329c';'2751326o';'2751855j';'2751849m';'2910557b' ;'2821044o';'3192141g';'2760704k'];
azimuth=[75.4897 -10             -3       -3         5.3          6.5   11.8         11.9      12.30      13          33           74           74         104];
magitude=[ 2.5 3.5          1.2        1.5        1.6           2.4   1.7        2.9        1.3        1.2       3.1         3.0           2.5         2.2];

ratall=zeros(13,nel,5);
ratlowathll=zeros(13,nel,3);
%'2960421s';%'2960421s';%'2960752d';%'2960421s'%'2960752d';%'2960421s';%'2960752d';%'3120406g4.p';%
for iii=1:1
filen=fall(iii,:);%'2960421s';%'2960752d';%'2960421s'%'2960752d';%'2960421s';%'2960752d';%'3120406g4.p';%

if magitude(iii)>2.4
    filename=[filen '1.p'];
else
    filename=[filen '4.p'];
end
[xtext pr]=load_bbdata([filename '01' ],[],dur);
start=pr.t0;

for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
  
    [xtext pr]=load_bbdata([filename '0' num2str(i)],[],dur);

    else
  
    [xtext pr]=load_bbdata([filename num2str(i)],[],dur);
   
    end
    
    tmp=pr.t0;
    if tmp(6)>start(6)
        start=tmp;
    end
 
end
k=1;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2740502k4
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);2821044o4.p
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x(k,:) pr]=load_bbdata([filename '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x(k,:) pr]=load_bbdata([filename num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x(k,:)=x(k,:)-mean(x(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end



if magitude(iii)>2.6
    for i=1:nel
        intx=zeros(1,length(x(1,:)));
        for j=1:length(x(i,:))
            intx(j)=sum(x(i,1:j));
        end
        x(i,:)=intx;
    end
end













          %[s c ph ci phi]=cmtm(x(1,:),x(2,:),0.005,1.5,0,0,1);




window=0.5;
ts1=1.5;
xx0=x(1,ts1*sr:(ts1+window)*sr);
lt=600;
tou=linspace(-0.5,0.5,lt);

figure(6);
for i=1:nel
subplot(nel,1,i);
plot(t,x(i,:),'r.');
xlim([1 3])
%xlim([ts1 ts1+window]);
end


for i=1:nel
    for j=1:lt
        clear tip xx;
        ts2=ts1+tou(j);
        tip=linspace(ts2,ts2+window,window*sr+1);
        xx=interp1(t,x(i,:),tip);
        xcr(i,j)=sum(xx0.*xx)/sqrt((sum(xx0.^2)*sum(xx.^2)));
    end
end
figure(7);
for i=1:nel
    subplot(12,1,i);
    plot(tou,xcr(i,:));
end
yy1=zeros(1,nel);
for i=1:nel
    for j=1:lt
        if xcr(i,j)==max(xcr(i,:))
            y1(i)=tou(j);
            mxcr(i)=xcr(i,j);
        end
    end
end
r=ones(nel,4);%coordinates of the array element
%r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
% r(:,1)=[-371 -298 -192 38 -23 0 23 177 213 260 98 -72 -131 -195]/1000;
% r(:,2)=[-302 -208 -284 -180 -13 0 -9 98 235 411 213 328 393 356]/1000;
% r(:,3)=[-25 -24 -16 -18 -4 0 2 18 12 3 -1 -16 -4 -2]/1000;
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
% y0=[4.545 4.502 4.491 4.434 4.421 4.413 4.333 4.310 4.20 4.325 4.42 4.41]
 %y1=[2.24 2.208 2.201 2.157 2.147 2.127 2.0935 2.0805 2.0465 2.081 2.198 2.183];%319 p
 % y1=[  2.2375    2.5984    2.2225    2.2024    2.3528    2.3227    2.9593    2.0370    2.3227    2.0621    2.0019    1.9768]

 %y1=[2.2205 2.205 2.1855 2.115 2.1054 2.145 2.085 2.0503 2.0155 2.0803 1.9753 1.966] ;%282 1
 %y1= [2.195 2.185 2.165  2.11 2.0855 2.1205 2.0552 2.0205 1.9705 2.0455 1.9303 1.9303] %282 2
% y1=[2.175 2.175 2.155 2.0904 2.095 2.121 2.0605 2.0145 1.9653 2.03515 1.9255 1.885]%291
% y2=[2.27 2.24 2.232 2.182 2.174 2.169 2.1215 2.108 2.077 2.1095 2.227 2.213];
% y3=[2.125 2.115 2.130 2.091 2.085 2.090 2.065 2.040 1.995 2.011 1.875 1.845];
% dt=[    0.0100    0.0116   -0.0090   -0.0081   -0.0064   -0.0087 -0.0032 -0.0059    0.0116    0.0319   -0.0333    0.0096] 319
% dt=[   -0.0187   -0.0213    0.0402   -0.0072   -0.0085   -0.0010 0.0145   -0.0075   -0.0179   -0.0008    0.0304   -0.0022] 282
%  dt=[-0.0152   -0.0227    0.0357   -0.0092   -0.0072    0.0002    0.0157  -0.0074   -0.0212    0.0028    0.0338   -0.0053] 291
%  dt=[   -0.0688   -0.0242    0.1026    0.0805    0.0444    0.0263   -0.2182   -0.0926    0.1066    0.1914   -0.3274    0.1795] 276
%  dt=[   -0.0214   -0.0219    0.0371   -0.0034   -0.0012   -0.0010  0.0133   -0.0110   -0.0220    0.0020    0.0333   -0.0040 ]274

% vp=0.9;
% %y1=[4.227 4.192 4.165 4.102 4.105 4.077 4.398 3.942 3.857 3.986 4.002 4.04];
% %y1=[4.232 4.187 4.165 4.108 4.092 4.076 3.992 3.942 3.857 3.988 4.002 4.036];
% xx=r\y1';
% 
% theta=atan2(xx(1),xx(2))*180/pi;
% phi=acos(xx(3)*vp)*180/pi;
% t0=xx(4);
% for i=1:nel
% y10(i)=xx(4)+xx(1:3)'*r(i,1:3)';
% end
% 
% dt=y10-y1
% figure(9);
% ys=sort(y1)
% for i=1:nel
%     for j=1:nel
%         if y1(j)==ys(i)
%             yss(i)=y10(j);
%             ydd(i)=y1(j);
%             if j>=4
%             id(i)=j+1;
%             else
%                 id(i)=j;
%             end
%         end
%     end
% end
% subplot(2,1,1);
% plot(ydd,yss,'-*');
% for i=1:nel
% text(ydd(i)+0.001,yss(i),num2str(id(i)));
% 
% end
% subplot(2,1,2);
% plot(dt,'-*');
% y10
% y1
% xx
% 1/sqrt(sum(xx(1:3).^2))
% 
% figure(7);
% for i=1:nel
% subplot(nel,1,i);
% ttt=linspace(1.5,2.5,sr+1);
% plot(ttt,x(i,round((1.5+y1(i))*sr):round((1.5+y1(i))*sr)+1*sr),'r.');
% xlim([1.5 2.5]);
% end












if magitude(iii)>2.4
    fileE=[filen '2.p'];
else
    fileE=[filen '5.p'];
end

k=1;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2740502k4
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);2821044o4.p
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [xe(k,:) pr]=load_bbdata([fileE '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [xe(k,:) pr]=load_bbdata([fileE num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    xe(k,:)=xe(k,:)-mean(xe(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end


if magitude(iii)>2.6
    for i=1:nel
        intx=zeros(1,length(xe(1,:)));
        for j=1:length(xe(i,:))
            intx(j)=sum(xe(i,1:j));
        end
        xe(i,:)=intx;
    end
end






if magitude(iii)>2.4
    fileN=[filen '3.p'];
else
    fileN=[filen '6.p'];
end

k=1;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2740502k4
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);2821044o4.p
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [xn(k,:) pr]=load_bbdata([fileN '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [xn(k,:) pr]=load_bbdata([fileN num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    xn(k,:)=xn(k,:)-mean(xn(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end


if magitude(iii)>2.6
    for i=1:nel
        intx=zeros(1,length(xn(1,:)));
        for j=1:length(xn(i,:))
            intx(j)=sum(xn(i,1:j));
        end
        xn(i,:)=intx;
    end
end





ttp1=1.5*sr;
ttp2=3*sr;

tts1=4*sr;
tts2=8*sr;
theta=azimuth(iii)-180;
xr=xn*cosd(theta)+xe*sind(theta);
xt=xn*sind(theta)+xe*cosd(theta);
% for page=1:4
% figure(17+page);
% oa=[1 2 3]+(page-1)*3;
% tlim1=1.6;
% tlim2=8;
% alim1=-1.2*1e3*max(xr(10,:));
% alim2=-alim1;
% for i=1:3
% subplot(9,1,3*i-2);
% plot(t,1e3*x(oa(i),:),'b');
% xlim([tlim1 tlim2]);
% ylim([alim1 alim2]);
% % rat(oa(i),1)=std(x(oa(i),ttp1:ttp2))/(std(xr(oa(i),ttp1:ttp2))^2+std(xt(oa(i),ttp1:ttp2))^2)^0.5;
% % ylabel(['v' num2str(rat(oa(i),1))]);
% subplot(9,1,3*i-1);
% plot(t,1e3*xr(oa(i),:),'r');
% % xlim([tlim1 tlim2]);
% % ylim([alim1 alim2]);
% rat(oa(i),2)=std(xr(oa(i),tts1:tts2))/std(xt(oa(i),tts1:tts2));
% ylabel(['r' num2str(rat(oa(i),2))]);
% subplot(9,1,3*i);
% plot(t,1e3*xt(oa(i),:),'g');
% xlim([tlim1 tlim2]);
% ylim([alim1 alim2]);
% % rat(oa(i),3)=std(xr(oa(i),ttp1:ttp2))/std(xt(oa(i),ttp1:ttp2));
% % ylabel(['t' num2str(rat(oa(i),3))]);
% 
% smin=1e8;
% smax=0;
% clear thmax thmin ss1 xtest;
% for th=0:360
%     xtest=xn(oa(i),:)*cosd(th)+xe(oa(i),:)*sind(th);
%     ss1=std(xtest(ttp1:ttp2));
%     if ss1>smax
%         smax=ss1;
%         thmax=th;
%     end
%     if ss1<smin
%         smin=ss1;
%         thmin=th;
%     end
% end
% 
% rat(oa(i),4)=thmax;
% rat(oa(i),5)=thmin;
% 
% 
% 
% end
% subplot(9,1,1);
% title([filen ' az=' num2str(round(theta)) 'sta' num2str(oa(1)) '-' num2str(oa(3)) ]);
% end
% 
% 

M=8;
    figure(17);
    title([filen ' vertical' ]);
    hold on;
    for i= 1:nel
    plot(t,x(i,:)/std(x(1,:))+i*M*2,'black');
    end
    xlim([1.5 inf]);

    figure(18);
    title([filen ' radial' ]);
    hold on;
    for i= 1:nel
    plot(t,xr(i,:)/std(xr(1,:))+i*M*2,'b');
    end
    xlim([1.5 inf]);
     figure(19);
     title([filen ' transverse' ]);
    hold on;
    for i= 1:nel
    plot(t,xt(i,:)/std(xt(1,:))+i*M*2,'r');
    end
    xlim([1.5 inf]);


% figure(22);
% plot(t,x(1,:)*4.5,'b');
% hold on;
% len=5;
% wid=0;
% hstf=ones(1,len);
% hstf(1:len-wid)=linspace(0,1,len-wid);
% stf=[hstf wrev(hstf)];
% xc=conv(x(1,:),stf);
% %plot(t,xc(len:length(xc)-len),'r');
% load 0752v;
% plot(t-0.2,aaa/1.2,'g');
% xlim([1.6 4]);




%     
% 
figure(26);

% D = fdesign.lowpass('Fp,Fst,Ap,Ast',3/100, 8/100, 1, 60);
%     butterbps = design(D,'equiripple');
%bps = fdesign.highpass(9.5, 10.8, 11.2,  .7, 7, 200);
%butterbps=butter(bps);
[B,A] = butter(4, 5/100,'high');

figure(27);
%[B,A]= butter(2, [12.5 14.5]/sr); %construct the filter for 8 to 12 hz
[H,F]=freqz(B,A,1024,sr);
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
yr=zeros(nel,dur*sr);
yt=zeros(nel,dur*sr);
yn=zeros(nel,dur*sr);
ye=zeros(nel,dur*sr);
for i=1:nel
    y(i,:)= filter(B,A,x(i,:));
    yr(i,:)= filter(B,A,xr(i,:));
    yt(i,:)= filter(B,A,xt(i,:));
    yn(i,:)= filter(B,A,xn(i,:));
    ye(i,:)= filter(B,A,xe(i,:));
end
y1 = filter(B,A,x(1,:));

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



M=8;
    figure(20);
    title([filen ' vertical' ]);
    hold on;
    for i= 1:nel
    plot(t,y(i,:)/std(y(1,:))+i*M,'black');
    end
    xlim([1.5 inf]);

    figure(21);
    title([filen ' radial' ]);
    hold on;
    for i= 1:nel
    plot(t,yr(i,:)/std(yr(1,:))+i*M,'b');
    end
    xlim([1.5 inf]);
     figure(22);
     title([filen ' transverse' ]);
    hold on;
    for i= 1:nel
    plot(t,yt(i,:)/std(yt(1,:))+i*M,'r');
    end
    xlim([1.5 inf]);








% ratlow=zeros(nel,3);
% for page=1:4
% figure(30+page);
% oa=[1 2 3]+(page-1)*3;
% tlim1=1.6;
% tlim2=8;
% alim1=-1.5*1e3*max(yr(10,:));
% alim2=-alim1;
% for i=1:3
% subplot(9,1,3*i-2);
% plot(t,1e3*y(oa(i),:),'b');
% xlim([tlim1 tlim2]);
% ylim([alim1 alim2]);
% ratlow(oa(i),1)=std(y(oa(i),ttp1:ttp2))/(std(yr(oa(i),ttp1:ttp2))^2+std(yt(oa(i),ttp1:ttp2))^2)^0.5;
% ylabel(['v' num2str(ratlow(oa(i),1))]);
% subplot(9,1,3*i-1);
% plot(t,1e3*yr(oa(i),:),'r');
% xlim([tlim1 tlim2]);
% ylim([alim1 alim2]);
% ratlow(oa(i),2)=std(yr(oa(i),tts1:tts2))/std(yt(oa(i),tts1:tts2));
% ylabel(['r' num2str(ratlow(oa(i),2))]);
% subplot(9,1,3*i);
% plot(t,1e3*yt(oa(i),:),'g');
% xlim([tlim1 tlim2]);
% ylim([alim1 alim2]);
% ratlow(oa(i),3)=std(yr(oa(i),ttp1:ttp2))/std(yt(oa(i),ttp1:ttp2));
% ylabel(['t' num2str(ratlow(oa(i),3))]);
% 
% smin=1e8;
% smax=0;
% clear thmax thmin ss1 ytest;
% for th=0:360
%     ytest=yn(oa(i),:)*cosd(th)+ye(oa(i),:)*sind(th);
%     ss1=std(ytest(ttp1:ttp2));
%     if ss1>smax
%         smax=ss1;
%         thmax=th;
%     end
%     if ss1<smin
%         smin=ss1;
%         thmin=th;
%     end
% end
% 
% ratlow(oa(i),4)=thmax;
% ratlow(oa(i),5)=thmin;
% 
% end
% subplot(9,1,1);
% title([filen ' az=' num2str(round(theta)) 'sta' num2str(oa(1)) '-' num2str(oa(3)) ]);
% end
% 
% 
% 
% ratall(iii,:,:)=rat;
% ratlowall(iii,:,:)=ratlow;
end


% 
 figure(50)
% subplot(3,1,3);
% plot(t,x(2,:));

subplot(3,1,1);
for i=1:length(x(2,:))
intx(i)=sum(x(2,1:i));
end
plot(t,intx);
subplot(3,1,2);
[B,A] = butter(4,1/100,'high');
plot(t,filter(B,A,intx));