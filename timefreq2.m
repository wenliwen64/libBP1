%close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy ];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
%r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=10;
t=linspace(0,dur,dur*sr);
x0=zeros(nel,dur*sr);
k=1;
filename='2821044o4.p';
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
start(6)=start(6)+300;

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

x1=x0;
figure(1);
plot(t,x2(1,:),t,x0(1,:),'r');


L=256/200;
ta1=1.5;
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
  sd=[   -0.0047   -0.0274   -0.0421   -0.0955   -0.1022   -0.1035   -0.1556   -0.1810   -0.2130   -0.1529   -0.1359   -0.1202]%282
  %sd=[   -0.0047   -0.0234   -0.0247   -0.0848   -0.0968   -0.0942   -0.1449   -0.1756   -0.2184   -0.1649   -0.1716   -0.1663]%291
for j=1:nel
    coe(j)=std(x1(1,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+0.04*sr-1))/std(x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+0.04*sr-1));
end
coe=[    1.0000    0.9413    0.6623    0.6512    0.7380    0.7968    0.8275    0.8153    0.7940    0.5561    0.7209    0.6309]
coe=ones(1,nel);
for j=1:nel
    xsum=xsum+coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    xorin(j,:)=coe(j)*x1(j,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+L*sr-1);
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
figure(9);
        w1= TFRBJ([xorin(1,:)' xorin(2,:)'],1:L*sr,200);
        ww=w1*0;
        tfRxx=zeros(sr,sr*L,nel,nel);
for i=1:nel
    for j=1:nel
        w= TFRBJ([xorin(i,:)' xorin(j,:)'],1:L*sr,200);
        for ii=1:sr
            for jj=1:sr*L
                tfRxx(ii,jj,i,j)=w(ii,jj);
            end
        end
    end
end
Rxx=zeros(nel,nel);
for ii=1:sr
    for jj=1:sr*L
        for i=1:nel
            for j=1:nel
                Rxx(i,j)=tfRxx(ii,jj,i,j);
            end
        end
        ww(ii,jj)=sum(eig(Rxx));
    end
end


%w1=tfrwv(x(1,:)',1:L,200);
%w2=tfrwv(x(12,:)',1:L,200);
%W=abs(w).^2./(abs(w1).*abs(w2));
%pcolor(W.^0.5);
pcolor(abs(ww));
colorbar;
ylim([0 100])
colorbar;
shading interp;

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


kkk=0;
%w=linspace(0,sr,window*sr);
for i=10:1:20
    kkk=kkk+1;
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    %Rxx=abs(X(:,i))*abs(X(:,i))';
       % Rxx=X(:,i)*X(:,i)';
   % R21a(i-fli+1)=Rxx(1,1);
%     Rxx=zeros(nel,nel);
    for j=1:nel
        for k=1:nel
        Rxx(j,k)=sum(tfRxx(i,100:150,j,k));
        end
    end
%     for j=1:nel
%         Rxx(j,j)=0.01;
%     end
    %R21b(i-fli+1)=Rxx(1,1);
    [Uv,A]=eig(Rxx);
    As=zeros(nel,nel);
    un=zeros(nel,nel);
    us=zeros(nel,nel);
    M=1;
    un(:,2:nel)=Uv(:,2:nel);
    us(:,nel:nel)=Uv(:,nel:nel);
    for iii=nel:nel
    As(iii,iii)=1/A(iii,iii);
    end
    Un=un*un';
    Us=us*us';
    %wi=fspec(i+fli-1);
        vi=i;
        %wi=s1(i);
        wi=1;%ww(kkk);
    for p=1:ps
        for q=1:qs
            
            a=v(r,kv(vi,ux(p),uy(q)));
            Pm(p,q)=Pm(p,q)+(wi*(a'*a)/(a'*Un*a));
        end
    end
end
% figure(13)
% subplot(2,1,1);
% plot(s1(fli:fhi),R21a,'b*');
% subplot(2,1,2);
% plot(s(fli:fhi),R21b,'b*');

figure(12);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
contourf(Xm,Ym,real(Pm),15);
%caxis([0 15]);
colorbar;
hold on;
plot(-0.186,-0.157,'white.');
plot(-0.06,-0.215,'white.');
hold off;
%title(['window start at' num2str(tts1) 'sec']);








