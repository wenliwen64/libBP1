close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
%r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=200;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2910557b4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+100;
window=100.5;
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
filename1='2821044o4.p';
[xtext pr]=load_bbdata([filename1 '13' ],[],dur);

start=pr.t0;
start(6)=start(6)+100;

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
coe291=[    1.0000    0.9241    0.6798    0.7913    0.8565    1.0021    1.3457    1.3643    1.4640    0.8563    0.8432    0.9188];%2910557b4
coe282=[    1.0000    0.8753    0.7198    0.7668    0.8069    0.8918    1.0063    0.6925    0.7285    0.5605    0.5551    0.5343];

x=zeros(nel,dur*sr);
%x=x1+x2;
for j=1:nel
   % x(j,:)=coe282(j)*x2(j,:)+coe291(j)*x1(j,:);
        x(j,:)=coe282(j)*x2(j,:)+coe291(j)*x1(j,:);
end

  ts1=0.4;
  uxx=0.3;
  uyy=0.2;
  sd=-[uxx uyy ]*r';
  
  
  
  ts2=0.4;
  uxx1=0.2;
  uyy1=0.2;
  sd1=-[uxx1 uyy1 ]*r';
  
  
 %ts2=0.4;
  uxx2=0.6;
  uyy2=0.3;
  sd2=-[uxx2 uyy2 ]*r';
  %sd=zeros(1,nel);
  tt=linspace(0,dur,sr*dur);
  ttt=linspace(0,window,window*sr);
  x1(1,:)=sin(2*pi*tt+0.947).*sin(2*pi*20*tt);
  %x2(1,:)=sin(2*pi*1.1*tt+3).*sin(2*pi*40*tt);
for j=1:nel
    %y(j,1:window*sr)=x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+2*x2(1,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
   y(j,1:window*sr)=x1(1,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr)+0*std(x1(1,:))*rand(1,window*sr);
end
%y=x;

figure(15);

plot(tt,x1(1,:),'r',tt,x2(1,:),'b');



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


for tl=2.2:0.1:2.2
clear X;

th=tl+2;
fl=20;
fh=20;
% s=linspace(0,sr,0.5*sr);
% for i=1:nel
%     X(i,:)=fft(y(i,tl*sr:th*sr-1));
%     
% end
 s1=linspace(0,sr,0.5*sr);
          [s c ph ci phi]=cmtm1(y(1,tl*sr:th*sr-1),y(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
     
         % fspec=linspace(0,sr/2,(th-tl)*sr/2);
          fli=round(interp1(s,1:length(s),fl));
          fhi=round(interp1(s,1:length(s),fh));
          ff=round(interp1(s,1:length(s),[3          13    23    ]));
          ww=                            [1          1      1      ];
          S=zeros(length(s),nel,nel);
for i=1:nel
    for j=1:nel
          [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1)-mean(y(i,tl*sr:th*sr-1)),y(j,tl*sr:th*sr-1)-mean(y(j,tl*sr:th*sr-1)),3);
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


kkk=0;
%w=linspace(0,sr,window*sr);
for i=fli:5:fhi
    kkk=kkk+1;
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    %Rxx=abs(X(:,i))*abs(X(:,i))';
     for k=1:20
        for j=1:nel
            X(j,:)=fft(y(j,tl*sr+(k-1)*0.5*sr:tl*sr+(k)*0.5*sr));
        end
         Rxx=Rxx+X(:,i)*X(:,i)';
    end
   % R21a(i-fli+1)=Rxx(1,1);
%     Rxx=zeros(nel,nel);
%     for j=1:nel
%         for k=1:nel
%         Rxx(j,k)=S(i,j,k);
%         end
%     end
%     %R21b(i-fli+1)=Rxx(1,1);
    [Uv,A]=eig(Rxx);
    As=zeros(nel,nel);
    un=zeros(nel,nel);
    us=zeros(nel,nel);
    M=5;
    un(:,1:M)=Uv(:,1:M);
    us(:,nel:nel)=Uv(:,nel:nel);
    for iii=nel:nel
    As(iii,iii)=1/A(iii,iii);
    end
    Un=un*un';
    Us=us*us';
    %wi=fspec(i+fli-1);
        vi=s(i);
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
end
