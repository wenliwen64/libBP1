close all;
clear all;
figure(1);
sr = 200; %sampling rate /sec

%t = linspace(0,1,sr);
%x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
dur=12;
nel=12;
x1=zeros(nel,dur*sr);
k=1;
for i=1:13
    if i==4
        continue;
    end
    if i<=9
    %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);
%     [x(k,:) pr]=load_bbdata(['3192141g6.p0' num2str(i)],[],dur);
    [x1(k,:) pr]=load_bbdata(['2910557b6.p0' num2str(i)],[],dur);
    else
    [x1(k,:) pr]=load_bbdata(['2910557b6.p' num2str(i)],[],dur);
    end
    x1(k,:)=x1(k,:)-mean(x1(k,:));
    k=k+1;
end
x2=zeros(nel,dur*sr);
k=1;
filename1='3192141g4.p';
[xtext pr]=load_bbdata([filename1 '01' ],[],dur);

start=pr.t0;
start(6)=start(6)+200;

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
x=x2;
%x=xx(1001:1500);
t=linspace(0,dur,dur*sr);
%subplot(3,2,1);                   
%[B,A]= cheby1(3,0.9, [16.9 23.7]/sr); %construct the filter for 8 to 12 hz
[B,A]= cheby1(3,0.9, [10 12]/sr); %construct the filter for 8 to 12 hz
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
for i=1:nel
    y(i,:)= filter(B,A,x(i,:));
end
y1 = filter(B,A,x(1,:));
y2= filter(B,A,x(2,:));
y3= filter(B,A,x(3,:));
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
figure(3);
for i=1:nel
subplot(nel,1,i)
plot(t,y(i,:),'b');
end
% subplot(3,1,2);
% plot(t,x(2,:),'b');
% subplot(3,1,3);1
% plot(t,x(2,:),'b');
%xlim([1 2]);
%subplot(3,2,5);                  
%zeros=roots(B);
%poles=roots(A);
%zplane(zeros,poles);
%xlabel('Real part');
%ylabel('imaginary');
%title('zeros and poles');
window=4.5;
step=5;
k=1;


time=window;
kkk=0;
for tts=0:0.25:9
    
%for ts1=0:step:dur-window

   clear yy;
    yy=zeros(nel,window*sr);
    yy=x(:,tts*sr+1:(tts+window)*sr);
    
    kkk=kkk+1;
M=3;%number of signals
k=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.3^2-ux^2-uy^2)^0.5];% wave vector
uv=@(ux,uy) [ux uy (1/0.3^2-ux^2-uy^2)^0.5];
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
r(:,2)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,1)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
 % sd=[ -0.0214   -0.0219    0.0371   -0.0034   -0.0012   -0.0010  0.0133 -0.0110   -0.0220    0.0020    0.0333   -0.0040];%274
  sd=[-0.0152   -0.0227    0.0357   -0.0092   -0.0072    0.0002    0.0157  -0.0074   -0.0212    0.0028    0.0338   -0.0053] ;%291

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




ps=30;
qs=30;
XX=zeros(ps,qs);
YY=zeros(ps,qs);

uxx=linspace(-1,1,ps);
uyy=linspace(-1,1,qs);
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

ts1=2;
wd=0.5;
for p=1:ps
    for q=1:qs
        clear uuv;
        uuv=uv(uxx(p),uyy(q));
        
        cc=zeros(nel,nel);
        for i=1:nel
            clear xxi;
             xxi=X(i,ts1*sr:(ts1+wd)*sr);
            for j=1:nel
                if i==j
                    continue;
                end
               
                clear xxj tip;
                touij=(rij(i,j,1)*uuv(1)+rij(i,j,2)*uuv(2)+rij(i,j,3)*uuv(3))+(sd(i)-sd(j));
               ts2=ts1+touij;
               tip=linspace(ts2,ts2+wd,wd*sr+1);
               xxj=interp1(t,X(j,:),tip);
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