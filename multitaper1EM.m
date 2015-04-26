%close all;
clear all;


figure(1);
%M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy ]%(1/2^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
%r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=12.5;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2821044o4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+100;
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
start(6)=start(6)+300;

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
%x=x1+x2;
%coe=[    0.9500    0.9413    0.6623    0.6512    0.7380    0.7968    0.8275    0.8153    0.7940    0.5561    0.7209    0.6309];
coe291=[    1.0000    0.9241    0.6798    0.7913    0.8565    1.0021    1.3457    1.3643    1.4640    0.8563    0.8432    0.9188];%2910557b4
coe282=[    1.0000    0.8753    0.7198    0.7668    0.8069    0.8918    1.0063    0.6925    0.7285    0.5605    0.5551    0.5343];
coe=ones(1,nel);
for j=1:nel
    x(j,:)=coe291(j)*x2(j,:);%+coe282(j)*x1(j,:);
end

  ts1=0.4;
  uxx=-0.25;
  uyy=-0.20;
  sd=zeros(1,nel);
  sd=-[uxx uyy ]*r';%(1/2^2-uxx^2-uyy^2)^0.5]*r';
  
  uxx1=-0.16;
  uyy1=-0.04;
  sd1=zeros(1,nel);
  sd1=-[uxx1 uyy1]*r';% (1/2^2-uxx1^2-uyy1^2)^0.5]*r';
  
 y=zeros(nel,window*sr);
 figure(14)
% for j=1:nel
%     y(j,1:window*sr)=2*x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1)+x2(10,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
%     subplot(nel,1,j);
%     plot(y(j,:));
% end
y=x;

figure(15);
tt=linspace(0,dur,sr*dur);
plot(tt,x1(9,:),'r',tt,x2(9,:),'b',tt,x3(9,:),'g');
%xlim([1.5 4.5]);

for tl=2.2:0.1:2.2
clear X;

th=tl+1;
fl=12;
fh=22;
%s=linspace(0,sr,0.5*sr);
figure(16)
for j=1:nel
   
    subplot(nel,1,j);
    plot(tt,x(j,:),tt,coe282(j)*x1(j,:),'r',tt,coe291(j)*x2(j,:),'g');
    xlim([tl th]);
    %ylim([-0.02 0.02]);
end


for i=1:nel
    X(i,:)=fft(y(i,tl*sr:th*sr-1)-mean(y(i,tl*sr:th*sr-1)));
    
end
 
          


[s c ph ci phi]=cmtm1(y(1,tl*sr:th*sr-1),y(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
fli=round(interp1(s,1:length(s),fl));
fhi=round(interp1(s,1:length(s),fh));


fran=1;



M=4;
Mw=[1/2 1/2];
N=60;
uran=0.2;

%uold=-[0.21 0.09 ;0.23 0.14];

A=zeros(length(s),nel,M);
A1=zeros(nel,M);
S=zeros(M,length(s));
Sn=zeros(M,length(s));
Xl=zeros(nel,length(s));
% uxo=zeros(1,M);
% uyo=zeros(1,M);
uxo=-rand(1,M)*0.3
uyo=-rand(1,M)*0.3

%              uxo(1)=-0.1697;
%              uyo(1)=-0.1179;
%              uxo(2)=-0.06;
%              uyo(2)=-0.21;
flag=0;
minuxn=zeros(1,M);
minuyn=zeros(1,M);
Mmax=zeros(1,M);
for iter=1:100
    uxn=zeros(1,M);
    uyn=zeros(1,M);
    for m=1:M
        clear Xl;
        for i=fli:fran:fhi
         a=v(r,kv(s(i),uxo(m),uyo(m)));
         A1(:,:)=A(i,:,:);
         Xl(:,i)=a*S(m,i)+1/M*(X(:,i)-A1*S(:,i));
        end
        
        uxs=linspace(uxo(m)-uran,uxo(m)+uran,N+1);
        uys=linspace(uyo(m)-uran,uyo(m)+uran,N+1);
        Mtestmax=0;
        for j=1:N
            for k=1:N
                Mtest=0;
                for i=fli:fran:fhi
                    a=v(r,kv(s(i),uxs(j),uys(k)));
                    Mtest=Mtest+a'*Xl(:,i)*Xl(:,i)'*a/(a'*a);
                end
                
                if Mtest>Mtestmax
                    uxn(m)=uxs(j);
                    uyn(m)=uys(k);
                    Mtestmax=Mtest;
                end
            end
        end
         Mmax(m)=Mtestmax;
%         
%               uxn(1)=-0.1655;
%              uyn(1)=-0.0979;
         
        for i=fli:fran:fhi
            a=v(r,kv(s(i),uxn(m),uyn(m)));
            Sn(m,i)=a'*Xl(:,i)/(a'*a);
        end
    end
    for m=1:M
        for i=fli:fran:fhi
        A(i,:,m)=v(r,kv(s(i),uxn(m),uyn(m)));
        end
    end
    S=Sn;
     if sum((uxo-uxn).^2)==0&&sum((uyo-uyn).^2)==0
         Mmax
        break
         %flag=flag+1;
    else
        flag=0;
    end
    if flag==1
       
        if sum((uxn-minuxn).^2)^0.5<0.02&&sum((uyn-minuyn).^2)^0.5<0.02
            Mmax
            break
        end
        minuxn=uxn;
        minuyn=uyn;
%         uxn=uxn+sign(rand(1,M)-0.5)*0.03;
%         uyn=uyn+sign(rand(1,M)-0.5)*0.03;
    end
    uxo=uxn;
    uyo=uyn;
    
    
uxo
uyo
end


figure(5);
hold on;
plot(uxo,uyo,'b-O');
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);
end