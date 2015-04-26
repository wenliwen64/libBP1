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
start(6)=start(6);

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

x1=x0+x2;
figure(1);
plot(t,x2(1,:),t,x0(1,:),'r');
for ii=1:1
L=3.0;
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
plot(tt,xres(j,:),'r',tt, xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b',tt,xsum,'g')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');

xlim([0 5.0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts1=0.1;
  uxx=-0.6;
  uyy=0.3;
 % sd=-[uxx uyy (1/0.7^2-uxx^2-uyy^2)^0.5]*r';
 sd=zeros(1,nel);
for j=1:nel
    %y(j,1:window*sr)=x(j,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
y=x;

% figure(15);
% tt=linspace(0,dur,sr*dur)
% plot(tt,x1(1,:),'r',tt,x2(1,:),'b');
% 


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


for tl=1.0:0.1:1.0
clear X;

th=tl+0.5;
fl=3;
fh=20;
% s=linspace(0,sr,0.5*sr);
for i=1:nel
    X(i,:)=fft(y(i,tl*sr:th*sr-1));
    
end
 s1=linspace(0,sr,0.5*sr);
          [s c ph ci phi]=cmtm1(y(1,tl*sr:th*sr-1),y(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
     
         % fspec=linspace(0,sr/2,(th-tl)*sr/2);
          fli=round(interp1(s,1:length(s),fl));
          fhi=round(interp1(s,1:length(s),fh));
          S=zeros(length(s),nel,nel);
for i=1:nel
    for j=1:nel
          [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1),y(j,tl*sr:th*sr-1),2);
          for k=1:length(s)
              S(k,i,j)=c(k);
          end
    end
end

ps=200;
qs=200;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs);
ux=linspace(-0.7,0.7,ps);
uy=linspace(-0.7,0.7,qs);

for q=1:qs
    Xm(:,q)=ux;
end
for p=1:ps
    Ym(p,:)=uy';
end



%w=linspace(0,sr,window*sr);
for i=fli:fhi
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
    %R21b(i-fli+1)=Rxx(1,1);
    [Uv,A]=eig(Rxx);
    
    un=zeros(nel,nel);
    un(:,1:nel-2)=Uv(:,1:nel-2);
    Un=un*un';
    %wi=fspec(i+fli-1);
        vi=s1(i);
        %wi=s1(i);
        wi=1;
    for p=1:ps
        for q=1:qs
            
            a=v(r,kv(vi,ux(p),uy(q)));
            Pm(p,q)=Pm(p,q)+wi*(a'*Un*a)/(a'*a);
        end
    end
end
% figure(13)
% subplot(2,1,1);
% plot(s1(fli:fhi),R21a,'b*');
% subplot(2,1,2);
% plot(s(fli:fhi),R21b,'b*');

figure(20);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
contourf(Xm,Ym,real(Pm),25);
%caxis([0 15]);
colorbar;
%title(['window start at' num2str(tts1) 'sec']);
figure(13);
clear minpm;
minpm=min(min(Pm));
for p=1:ps
    for q=1:qs
        if Pm(p,q)==minpm
            pm(ii)=ux(p);
            qm(ii)=uy(q);
            break;
        end
        
    end
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
figure(17)
plot(pm,qm,'.');
%xlim([-0.3 0.0]);
%ylim([-0.3 0.0]);