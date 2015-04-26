close all;
clear all;


figure(1);
%M=1;
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
x=x1+2*x2;


  ts1=0.2;
  uxx=-0.04;
  uyy=-0.04;
  sd=zeros(1,nel);
  sd=-[uxx uyy (1/2^2-uxx^2-uyy^2)^0.5]*r';
  
  uxx1=-0.16;
  uyy1=-0.16;
  sd1=zeros(1,nel);
  sd1=-[uxx1 uyy1 (1/2^2-uxx1^2-uyy1^2)^0.5]*r';
  
%  y=zeros(nel,window*sr);
%  figure(14)
% for j=1:nel
%     y(j,1:window*sr)=x(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1)+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
%     subplot(nel,1,j);
%     plot(y(j,:));
% end
y=x;

figure(15);
tt=linspace(0,dur,sr*dur);
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


for tl=1.5:0.1:1.5
clear X;

th=tl+3;
fl=10;
fh=10;
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
          [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1),y(j,tl*sr:th*sr-1),4);
          for k=1:length(s)
              S(k,i,j)=c(k);
          end
    end
end

ps=20;
qs=20;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs,ps,qs);
ux=linspace(-0.24,-0.20,ps);
uy=linspace(-0.13,-0.06,qs);


% t(-0.2172, -0.1230,'white.');
% plot(-0.06,-0.215,'white.');
% plot(-0.2188,-0.0927,'white.');



for q=1:qs
    Xm(:,q)=ux;
end
for p=1:ps
    Ym(p,:)=uy';
end


%     %R21b(i-fli+1)=Rxx(1,1);
%     [Uv,A]=eig(Rxx);
%     As=zeros(nel,nel);
%     un=zeros(nel,nel);
%     us=zeros(nel,nel);
%     M=5;
%     un(:,1:nel-7)=Uv(:,1:nel-7);
%     us(:,nel:nel)=Uv(:,nel:nel);
%     for iii=nel:nel
%     As(iii,iii)=1/A(iii,iii);
%     end
%     Un=un*un';
%     Us=us*us';
%     %wi=fspec(i+fli-1);
%         vi=s(i);
%         %wi=s1(i);
%         wi=1;%ww(kkk);
%     for p=1:ps
%         for q=1:qs
%             
%             a=v(r,kv(vi,ux(p),uy(q)));
%             Pm(p,q)=Pm(p,q)+(wi*(a'*Us*a)/(a'*Un*a));
%         end
%     end
M=2;
N=2;
Rxx=zeros(nel,nel);
unew=zeros(M,N);
uold=zeros(M,N);
uold=-[0.21 0.09 ;0.23 0.14];
utest=zeros(M,N);
A=zeros(length(s),nel,M);
A1=zeros(nel,M);
pmtesti=zeros(1,length(s));
for iter=1:1000
    for m=1:M
        for n=1:N
            utest=uold;
          
            pmtestold=0;
            uutest=linspace(uold(m,n)-0.02,uold(m,n)+0.02,21);
            for k=1:21
                pmtesti=zeros(1,length(s));
                utest(m,n)=uutest(k);
                if utest(1,1)==utest(2,1)||utest(1,2)==utest(2,2)
                    continue;
                end
                clear PA Ap PAp Rxx;
                for i=fli:10:fhi
                     A1=zeros(nel,M);
                    A(i,:,m)=v(r,kv(s(i),utest(m,1),utest(m,2)));
                    for j=1:nel
                        for jj=1:M
                            A1(j,jj)=A(i,j,jj);
                        end
                    end
                    Ap=(A1'*A1)\A1';
                    PA=A1*Ap;
                    PAp=eye(nel)-PA;
                    sigma=1/(nel-M)*trace(PAp);
                    
                    for j=1:nel
                        for jj=1:nel
                            Rxx(j,jj)=S(i,j,jj);
                        end
                    end
                    Psml=Ap*(Rxx-sigma*eye(nel))*Ap';
                    pmtesti(i)=trace(A1*Psml*A1'+sigma*eye(nel));
                    
                end
                pmtestnew=sum(pmtesti);
                if pmtestnew>pmtestold
                    pmtestold=pmtestnew;
                    unew(m,n)=utest(m,n);
                end
            end
        end
        
    end
    uold=unew;
     for i=fli:10:fhi
         for m=1:M
              A(i,:,m)=v(r,kv(s(i),uold(m,1),uold(m,2)));
         end
     end
end
uold
end