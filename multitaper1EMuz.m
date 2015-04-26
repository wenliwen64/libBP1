%close all;
clear all;


figure(1);
%M=1;
nel=12;
kv=@(f0,ux,uy,uz) 2*pi*f0*[ux uy uz];% wave vector
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
start(6)=start(6)+100;

for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2821044o4.p
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
x=x1+x2;


%   ts1=0.2;
%   uxx=-0.04;
%   uyy=-0.10;
%   sd=zeros(1,nel);
%   sd=-[uxx uyy (1/2^2-uxx^2-uyy^2)^0.5]*r';
%   
%   uxx1=-0.16;
%   uyy1=-0.10;
%   sd1=zeros(1,nel);
%   sd1=-[uxx1 uyy1 (1/2^2-uxx1^2-uyy1^2)^0.5]*r';
%   
%  y=zeros(nel,window*sr);
%  figure(14)
% for j=1:nel
%     y(j,1:window*sr)=1.2*x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1)+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
%     subplot(nel,1,j);
%     plot(y(j,:));
% end
y=x;

figure(15);
tt=linspace(0,dur,sr*dur);
plot(tt,x1(1,:),'r',tt,x2(1,:),'b',tt,x3(1,:),'g');
xlim([1.5 4.5]);

for tl=2:0.1:2
clear X;

th=tl+1;
fl=10;
fh=30;
%s=linspace(0,sr,0.5*sr);
for i=1:nel
    X(i,:)=fft(y(i,tl*sr:th*sr-1));
    
end
 
          


[s c ph ci phi]=cmtm1(y(1,tl*sr:th*sr-1),y(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
fli=round(interp1(s,1:length(s),fl));
fhi=round(interp1(s,1:length(s),fh));
fran=2;



M=2;
N=20;
uran=0.1;

%uold=-[0.21 0.09 ;0.23 0.14];

A=zeros(length(s),nel,M);
A1=zeros(nel,M);
S=zeros(M,length(s));
Sn=zeros(M,length(s));
Xl=zeros(nel,length(s));
% uxo=zeros(1,M);
% uyo=zeros(1,M);
uxo=-rand(1,M)*0.2;
uyo=-rand(1,M)*0.2;
uzo=rand(1,3);
for iter=1:20
    uxn=zeros(1,M);
    uyn=zeros(1,M);
    uzn=zeros(1,M);
    for m=1:M
        clear Xl;
        for i=fli:fran:fhi
         a=v(r,kv(s(i),uxo(m),uyo(m),uzo(m)));
         A1(:,:)=A(i,:,:);
         Xl(:,i)=a*S(m,i)+1/M*(X(:,i)-A1*S(:,i));
        end
        
        uxs=linspace(uxo(m)-uran,uxo(m)+uran,N+1);
        uys=linspace(uyo(m)-uran,uyo(m)+uran,N+1);
        uzs=linspace(uzo(m)-uran,uzo(m)+uran,N+1);
        Mtestmax=0;
        for j=1:N
            for k=1:N
                for kk=1:N
                    Mtest=0;
                    for i=fli:fran:fhi
                        a=v(r,kv(s(i),uxs(j),uys(k),uzs(kk)));
                        Mtest=Mtest+a'*Xl(:,i)*Xl(:,i)'*a/(a'*a);
                    end
                    
                    if Mtest>Mtestmax
                        uxn(m)=uxs(j);
                        uyn(m)=uys(k);
                        uzn(m)=uzs(kk);
                        Mtestmax=Mtest;
                    end
                end
            end
        end
        for i=fli:fran:fhi
            a=v(r,kv(s(i),uxn(m),uyn(m),uzn(m)));
            Sn(m,i)=a'*Xl(:,i)/(a'*a);
        end
    end
    for m=1:M
        for i=fli:fran:fhi
        A(i,:,m)=v(r,kv(s(i),uxn(m),uyn(m),uzn(m)));
        end
    end
    S=Sn;
    if sum((uxo-uxn).^2)==0&&sum((uyo-uyn).^2)==0
        break;
    end
    uxo=uxn;
    uyo=uyn;
    uzo=uzn;
    
uxo
uyo
uzo
    
end

figure(5);
hold on;
plot(uxo,uyo,'b-o');
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);
end