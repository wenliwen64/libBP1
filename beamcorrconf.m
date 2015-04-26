close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/2.0^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=10;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2821044o4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6);
window=2.5;
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
for ii=1:30
L=1.0;
ta1=2.05;
ta2=ta1+L;
xsum=zeros(1,L*sr);
yr=zeros(nel,L*sr);
xrr=zeros(nel,L*sr);

tt=linspace(ta1,ta2,L*sr);
sd=[   -0.0033   -0.0260   -0.0421   -0.0968   -0.1008   -0.1022   -0.1543   -0.1823   -0.2117   -0.1556   -0.1356   -0.1262];
for j=1:nel
    xsum=xsum+x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
xsum=xsum/nel;
xres=zeros(1,L*sr);
xr=x1;
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
    xxr(j,:)=abs(ifft(abs(yr(j,:)).*(cos(ph)+1i*sin(ph))));
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
plot(tt,xres(j,:),'r',tt, xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b',tt,xsum,'g')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');

xlim([0 3.0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window=4.5;
step=5;
k=1;


time=window;
kkk=0;
for tts=0.5:0.1:0.5
    
%for ts1=0:step:dur-window

   clear yy;
    yy=zeros(nel,window*sr);
    yy=x(:,tts*sr+1:(tts+window)*sr);
    
    kkk=kkk+1;
M=3;%number of signals
k=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.3^2-ux^2-uy^2)^0.5];% wave vector
uv=@(ux,uy) [ux uy (1/1.7^2-ux^2-uy^2)^0.5];
v=@(r,k) exp(1i*(r*k'));%steering vector
delay= @(x,y,the,c) (x/cos(the)+(y-tan(the)*x)*sin(the))/c;
d=0.05;%km
nel=12;
%km
r=zeros(nel,3);%coordinates of the array element
%r(:,1)=([1:nel]*d-d)-(nel*d-d)/2;
% r(:,1)=[-371 -298 -192 38 -23 0 23 177 213 260 98 -72 -131 -195]/1000;
% r(:,2)=[-302 -208 -284 -180 -13 0 -9 98 235 411 213 328 393 356]/1000;
% r(:,3)=[-25 -24 -16 -18 -4 0 2 18 12 3 -1 -16 -4 -2]/1000;
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
  sd274=[ -0.0214   -0.0219    0.0371   -0.0034   -0.0012   -0.0010  0.0133 -0.0110   -0.0220    0.0020    0.0333   -0.0040];%274
  sd291=[-0.0152   -0.0227    0.0357   -0.0092   -0.0072    0.0002    0.0157  -0.0074   -0.0212    0.0028    0.0338   -0.0053] ;%291

 %sd=[     0.0035    0.0000    0.0039   -0.0079   -0.0032   -0.0074    0.0057    0.0039    0.0010   -0.0003    0.0013   -0.0006];% 282
 sd=zeros(1,12);
  %sd=(sd274+sd291+sd282)/3;
 %time=2;%total time in seconds

samplerate=200;
t=linspace(0,time,time*samplerate);
L=time*samplerate;
fc=6.5;
rij=zeros(nel,nel,3);
for i=1:nel
    for j=1:nel
        for k=1:3
        rij(i,j,k)=r(i,k)-r(j,k);
        end
    end
end

clear X;
X=yy;




ps=100;
qs=100;
XX=zeros(ps,qs);
YY=zeros(ps,qs);

uxx=linspace(-0.3,-0.10,ps);
uyy=linspace(-0.2,-0.0,qs);
for q=1:qs
    XX(:,q)=uxx;
end
for p=1:ps
    YY(p,:)=uyy';
end
acc=zeros(ps,qs);
% 
% window=0.5;
% ts1=2;
% xx0=x(1,ts1*sr:(ts1+window)*sr);
% lt=300;
% tou=linspace(-0.4,0.4,lt);
% for i=1:nel
%     for j=1:lt
%         clear tip xx;
%         ts2=ts1+tou(j);
%         tip=linspace(ts2,ts2+window,window*sr+1);
%         xx=interp1(t,x(i,:),tip);
%         xcr(i,j)=sum(xx0.*xx)/sqrt((sum(xx0.^2)*sum(xx.^2)));
%     end
% end

ts1=1.2;
wd=0.5;
for p=1:ps
    for q=1:qs
        clear uuv;
        uuv=uv(uxx(p),uyy(q));
        if imag(uuv(3))~=0;
            continue;
        end
        cc=zeros(nel,nel);
        for i=1:nel
            clear xxi;
             xxi=X(i,ts1*sr:(ts1+wd)*sr);
            for j=1:nel
                if i==j
                    continue;
                end
               
                %clear xxj tip;
                touij=-(rij(i,j,1)*uuv(1)+rij(i,j,2)*uuv(2)+rij(i,j,3)*uuv(3))+(sd(i)-sd(j));
               ts2=ts1+touij;
%                tip=linspace(ts2,ts2+wd,wd*sr+1);
%                 xxj=interp1(t,X(j,:),tip);
                xxj=X(j,round(ts2*sr):round((ts2)*sr)+wd*sr);
                cc(i,j)=sum(xxi.*xxj)/sqrt((sum(xxi.^2)*sum(xxj.^2)));
                %cc(i,j)=sum(X(i,ts1:ts2).*X(j,ts1+touij:ts2+touij))/sqrt((sum(X(i,ts1:ts2).^2)*sum(X(j,ts1+touij:ts2+touij).^2)));
                
            end
        end
        acc(p,q)=sum(sum(cc))/(nel*nel-nel);
    end
end
figure(kkk+10);
contourf(XX,YY,acc,15);
colorbar;
end
maxacc=max(max(acc));
for p=1:ps
    for q=1:qs
        if acc(p,q)==maxacc
            pm(ii)=uxx(p);
            qm(ii)=uyy(q);
            break;
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
figure(17)
plot(pm,qm,'.');
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);