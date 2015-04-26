close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/2^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=20;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2821044o4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6);
window=5.5;
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
filename1='2910557b4.p';
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

x3=zeros(nel,dur*sr);
k=1;
filename1='3192141g4.p';
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
    [x3(k,:) pr]=load_bbdata([filename1 '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x3(k,:) pr]=load_bbdata([filename1 num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x3(k,:)=x3(k,:)-mean(x3(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end
x=zeros(nel,dur*sr);
x=x1;


  ts1=0.1;
  uxx=0.1;
  uyy=0.1;
  sd=-[uxx uyy (1/2^2-uxx^2-uyy^2)^0.5]*r';
 sd=zeros(1,nel);
for j=1:nel
    y(j,1:window*sr)=x(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
%y=x;

figure(15);
tt=linspace(0,dur,sr*dur);
plot(tt,x1(1,:),'r',tt,x2(1,:),'b',tt,2*x3(1,:),'g');

%  H=zeros(nel,M);
% for i=1:M
%     H(:,i)=v(r,k(fc,ux(i),uy(i)));
% end
% X=fft(real(H*S+U));

for tl=1.5:0.1:1.5
clear X;

th=tl+2;
fl=7;
fh=15;
% s=linspace(0,sr,0.5*sr);
% for i=1:nel
% X=zeros(nel,(th-tl)*sr);
% for i=1:nel
%      X(i,:)=fft(y(i,tl*sr:th*sr-1));
%      end
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
          [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1),y(j,tl*sr:th*sr-1),4);
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

 M=1;
 fr=10;
 Ar=zeros(nel,M);
 %umx=[-0.186];
 %umy=[-0.06];
%  umx=[-0.06];%[-0.186 -0.06];
%  umy=[-0.215];%[-0.157 -0.215];
umx=[0.1 -0.2188 -0.06 ];
umy=[ 0.1  -0.0927 -0.215];
 for i=1:M
     Ar(:,i)=v(r,kv(fr,umx(i),umy(i)));
 end
  Rgen=zeros(nel,nel);
for i=fli:1:fhi
    kkk=kkk+1;
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    %Rxx=abs(X(:,i))*abs(X(:,i))';
        %Rxx=X(:,i)*X(:,i)';
   % R21a(i-fli+1)=Rxx(1,1);
%     Rxx=zeros(nel,nel);
fi=s(i);
clear Ai Ui Ci Vi Ti;
Ai=zeros(nel,M);
for j=1:M
     Ai(:,j)=v(r,kv(fi,umx(j),umy(j)));
end

[Ui Ci Vi]=svd(Ar*Ai') ;
Ti=Vi*Ui';
    for j=1:nel
        for k=1:nel
        Rxx(j,k)=S(i,j,k);
        end
    end
  Rgen=Rgen+Ti*Rxx*Ti';  
end

    %R21b(i-fli+1)=Rxx(1,1);
    [Uv,A]=eigs(Rgen,nel,'LM');
    As=zeros(nel,nel);
    un=zeros(nel,nel);
    us=zeros(nel,nel);
    N=7;
    un(:,N:nel)=Uv(:,N:nel);
    us(:,1:2)=Uv(:,1:2);
    for iii=nel:nel
    As(iii,iii)=1/A(iii,iii);
    end
    Un=un*un';
    Us=us*us';
    %wi=fspec(i+fli-1);
        
        %wi=s1(i);
        wi=1;%ww(kkk);
    for p=1:ps
        for q=1:qs
            
            a=v(r,kv(fr,ux(p),uy(q)));
            Pm(p,q)=(wi*(a'*Us*a)/(a'*Un*a))^1;
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

   
  
plot(-0.2172, -0.1230,'white.');
plot(-0.06,-0.215,'white.');
plot(-0.2188,-0.0927,'white.');
hold off;
%title(['window start at' num2str(tts1) 'sec']);
end
