close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.70^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=10;
t=linspace(0,dur,dur*sr);
x0=zeros(nel,dur*sr);
k=1;
filename='2721715j1.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+100;
window=2.5;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2821044o4.p
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x0(k,:) pr]=load_bbdata([filename '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x0(k,:) pr]=load_bbdata([filename num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x0(k,:)=x0(k,:)-mean(x0(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end
x2=zeros(nel,dur*sr);
k=1;
filename1='2910557b4.p';
[xtext pr]=load_bbdata([filename1 '13' ],[],dur);

start=pr.t0;
start(6)=start(6)+100;

for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);filename1='2740502k4.p';
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

x1=1*x0;
figure(1);
plot(t,x2(1,:),t,x0(1,:),'r');


L=5.0;
ta1=0.8;
ta2=ta1+L;
xsum=zeros(1,L*sr);
yr=zeros(nel,L*sr);
xrr=zeros(nel,L*sr);
xorin=zeros(nel,L*sr);
tt=linspace(ta1,ta2,L*sr);
%sd=[   -0.0033   -0.0260   -0.0421   -0.0968   -0.1008   -0.1022   -0.1543   -0.1823   -0.2117   -0.1556   -0.1356   -0.1262];
 %sd=  [-0.0030   -0.0270   -0.0471   -0.0952   -0.1012   -0.1032   -0.1553   -0.1813   -0.2154   -0.1533   -0.1392   -0.1312]
 %sd=[     -0.0033   -0.0234   -0.0234   -0.0848   -0.0928   -0.0888   -0.1422   -0.1730   -0.2170   -0.1609   -0.1716   -0.1663]
  %sd=[     -0.0047   -0.0220   -0.0354   -0.0875   -0.0955   -0.0955   -0.1449   -0.1649   -0.1917   -0.1409   -0.1195   -0.1048] % 319
  %sd=[   -0.0047   -0.0274   -0.0421   -0.0955   -0.1022   -0.1035   -0.1556   -0.1810   -0.2130   -0.1529   -0.1359   -0.1202]%282
 % sd=[   -0.0047   -0.0234   -0.0247   -0.0848   -0.0968   -0.0942   -0.1449   -0.1756   -0.2184   -0.1649   -0.1716   -0.1663]%291
   sd=[ 0.0068   -0.0142   -0.0324   -0.0767   -0.0804   -0.0845   -0.1164   -0.1342   -0.1581   -0.1132   -0.0850   -0.0652];%272
for j=1:nel
    coe(j)=std(x1(1,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+0.04*sr-1))/std(x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+0.04*sr-1));
end
coe=[    1.0000    0.9413    0.6623    0.6512    0.7380    0.7968    0.8275    0.8153    0.7940    0.5561    0.7209    0.6309]
for j=1:nel
    xsum=xsum+coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    xorin(j,:)=x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    xr(j,:)=coe(j)*x1(j,:);
end
xsum=xsum/nel;
xres=zeros(1,L*sr);
%xr=x1;
for j=1:nel
    xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)=xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)-xsum;
end
for j=1:nel
    xres(j,:)=xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    yr(j,:)=fft(xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1));
end
for j=1:nel
    ph=2*pi*rand(1,L*sr);
    xxr(j,:)=ifft(abs(yr(j,:)).*(cos(ph)+1i*sin(ph)));
end
for j=1:nel
    xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)=xsum+xxr(j,:);
end
x=xr;

figure(7);
%plot(tt,xsum);
for j=1:nel
subplot(nel,1,j);
%ttt=linspace(1.5,2.5,sr+1);
%plot(tt,xres(j,:),'r',tt, coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b',tt,xsum,'g')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');
plot(tt,xorin(j,:));
%xlim([1.7 3.3]);
%ylim([-0.02 0.02]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ts1=0.1;
%   uxx=-0.6;
%   uyy=0.3;
%  % sd=-[uxx uyy (1/0.7^2-uxx^2-uyy^2)^0.5]*r';
%  sd=zeros(1,nel);
% for j=1:nel
%     %y(j,1:window*sr)=x(j,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
% end
% y=x;
%      Rxx=zeros(nel,nel);
% for i=1:length(tt)
%     Rxx=Rxx+xres(:,i)*xres(:,i)';
% end
% for i=1:nel
%     for j=1:nel
%         Rxx1(i,j)=Rxx(i,j)/sqrt(sum(xres(i,:).^2)*sum(xres(j,:).^2));
%     end
% end
% % figure(15);
% tt=linspace(0,dur,sr*dur)
% plot(tt,x1(1,:),'r',tt,x2(1,:),'b');
% 
figure(26);
bps = fdesign.bandpass(9.5, 10.8, 11.2, 12.5, 10, .7, 7, 200);
butterbps=butter(bps);
figure(27);
%[B,A]= butter(2, [12.5 14.5]/sr); %construct the filter for 8 to 12 hz
[B,A]= butter(3, 140/sr,'high'); %construct the filter for 8 to 12 hz
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
yorin=zeros(nel,L*sr);
for i=1:nel
    yorin(i,:)= filter(B,A,xorin(i,:));
end
figure(8);
for j=1:nel
subplot(nel,1,j);
%ttt=linspace(1.5,2.5,sr+1);
%plot(tt,xres(j,:),'r',tt, coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b',tt,xsum,'g')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');
plot(tt,xorin(j,:));
%xlim([1.7 3.3]);
%ylim([-0.02 0.02]);
end
% figure(12);
% %plot(t,x1(1,:),'b',t,x2(1,:),'r');
% Xtw=zeros(sr*dur,sr*dur);
% Yfw=zeros(sr*dur,sr*dur);
% for i=1:dur*sr
% Xtw(:,i)=linspace(0,dur,dur*sr);
% end
% for i=1:dur*sr
% Yfw(i,:)=linspace(-sr/2,sr/2,dur*sr);
% end
% w=mywigner(x(1,:)');
% pcolor(Yfw,Xtw,abs(w));
% shading interp;
% xlim([0 20]);
% 
% 
mm=70
Pm=zeros(25,mm);
kkk=0;
for tl=0.2:0.05:0.05*(mm-1)+0.2
clear X;
kkk=kkk+1
th=tl+0.25;
if tl==1.8
    th=2.0525;
end
fl=3;
fh=20;
% s=linspace(0,sr,0.5*sr);
% for i=1:nel
%     X(i,:)=fft(y(i,tl*sr:th*sr-1));
%     
% end
 s1=linspace(0,sr,0.5*sr);
          [s c ph ci phi]=cmtm1(yorin(1,tl*sr:th*sr-1),yorin(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
     
         % fspec=linspace(0,sr/2,(th-tl)*sr/2);
          fli=round(interp1(s,1:length(s),fl));
          fhi=round(interp1(s,1:length(s),fh));
          S=zeros(length(s),nel,nel);
for i=1:nel
    for j=1:nel
          [c ph ci ]=cmtm3(xorin(i,tl*sr:th*sr-1),xorin(j,tl*sr:th*sr-1),2);
          for k=1:length(s)
              S(k,i,j)=c(k);
          end
    end
end

% ps=200;
% qs=200;
% Xm=zeros(ps,qs);
% Ym=zeros(ps,qs);
% Pm=zeros(ps,qs);
% ux=linspace(-0.7,0.7,ps);
% uy=linspace(-0.7,0.7,qs);
% 
% for q=1:qs
%     Xm(:,q)=ux;
% end
% for p=1:ps
%     Ym(p,:)=uy';
% end

Coh=zeros(1,length(s));

%w=linspace(0,sr,window*sr);
for i=1:length(s)
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    %Rxx=abs(X(:,i))*abs(X(:,i))';
       % Rxx=X(:,i)*X(:,i)';
   % R21a(i-fli+1)=Rxx(1,1);
%     Rxx=zeros(nel,nel);
    for j=1:nel
        for k=1:nel
        Rxx(j,k)=S(i,j,k);
        end
    end
    Coh(i)=(sum(sum(Rxx))-trace(Rxx))/(nel*nel-nel);
end
 
%     figure(13)
%     subplot(mm,2,2*kkk-1);
%     plot(s,Coh);
%     ylim([0.3 1]);
%     subplot(mm,2,2*kkk);
%     ttt=linspace(tl,th,0.25*sr);
%     plot(ttt,xorin(1,tl*sr:th*sr-1));
%      xlim([tl th])
     Pm(:,kkk)=Coh;
     
%     %R21b(i-fli+1)=Rxx(1,1);
%     [Uv,A]=eig(Rxx);
%     
%     un=zeros(nel,nel);
%     un(:,1:nel-2)=Uv(:,1:nel-2);
%     Un=un*un';
%     %wi=fspec(i+fli-1);
%         vi=s1(i);
%         %wi=s1(i);
%         wi=1;
%     for p=1:ps
%         for q=1:qs
%             
%             a=v(r,kv(vi,ux(p),uy(q)));
%             Pm(p,q)=Pm(p,q)+wi*(a'*Un*a)/(a'*a);
%         end
%     end
% end
% % figure(13)
% % subplot(2,1,1);
% % plot(s1(fli:fhi),R21a,'b*');
% % subplot(2,1,2);
% % plot(s(fli:fhi),R21b,'b*');
% 
% figure(20);
% %contourf(X,Y,real(Pm),0.5:0.1:3.5);
% %caxis([0.5 3.5]);
% contourf(Xm,Ym,real(Pm),25);
% %caxis([0 15]);
% colorbar;
% %title(['window start at' num2str(tts1) 'sec']);
% figure(13);
% clear minpm;
% minpm=min(min(Pm));
% for p=1:ps
%     for q=1:qs
%         if Pm(p,q)==minpm
%             pm(ii)=ux(p);
%             qm(ii)=uy(q);
%             break;
%         end
%         
%     end
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% for i=1:25
%     for j=1:mm
%         if Pm(i,j)<0.5
%             Pm(i,j)=NaN;
%         end
%     end
% end

ux=linspace(1,1+0.05*(mm-1),mm);
uy=linspace(0,100,25);

for q=1:25
    Xm(:,q)=ux;
end
for p=1:mm
    Ym(p,:)=uy';
end

figure(17);
pcolor(Xm,Ym,Pm');
caxis([0 1]);
colorbar;
shading flat;
ylim([0 40]);
xlabel('time(s)');
ylabel('frequency(Hz)');
title('time-frequency coherece coefficient (0-1)');
%shading interp;
% figure(17)
% plot(pm,qm,'.');
% end
%xlim([-0.3 0.0]);
%ylim([-0.3 0.0]);