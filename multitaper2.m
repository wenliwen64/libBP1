%close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/1.7^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=8;
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

x2=zeros(nel,dur*sr);
k=1;
filename1='2821044o4.p';
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


% x=zeros(nel,dur*sr);
% x=x2;


  ts1=0.1;
  uxx=-0.6;
  uyy=0.3;
 % sd=-[uxx uyy (1/0.7^2-uxx^2-uyy^2)^0.5]*r';
 sd=zeros(1,nel);
for j=1:nel
    %y(j,1:window*sr)=x(j,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
%y=x;

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

% figure(2);
% %[TFR,T,F]=tfrspwv(x(1,:)',1:0.1:dur-1,40);\
% 
% w=tfrspwv([x(10,:)' x(11,:)'],1:L,200);
% w1=tfrspwv(x(1,:)',1:L,200);
% w2=tfrspwv(x(10,:)',1:L,200);
% W=abs(w)/(abs(w1).*abs(w2));
% pcolor(abs(w));
% colorbar;
% ylim([0 100])
% colorbar;
% shading interp;




clear X;
L=256;
N=256;
tl=0.40;
th=1.0;
fl=8;
fh=19;
fli=round(fl*N/sr);
fhi=round(fh*N/sr);
% s=linspace(0,sr,0.5*sr);
x=zeros(nel,L);
x=x2(:,1+4.5*sr:L+4.5*sr);
figure(3);
 w=tfrspwv([x(1,:)' x(11,:)'],1:L,N);
 pcolor(abs(w));
          
          S=zeros(nel,nel,N,L);
for i=1:nel
    for j=1:nel
          S(i,j,:,:)=tfrspwv([x(i,:)' x(j,:)'],1:L,N);
          
    end
end

ps=100;
qs=100;
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
     Rxx=zeros(nel,nel);
     for j=1:nel
         for k=1:nel
             for kk=tl*sr:th*sr-1
                 Rxx(j,k)=Rxx(j,k)+S(j,k,i,kk);
             end
         end
     end
    %R21b(i-fli+1)=Rxx(1,1);
    [Uv,A]=eig(Rxx);
    
    un=zeros(nel,nel);
    un(:,1:nel-1)=Uv(:,1:nel-1);
    Un=un*un';
    %wi=fspec(i+fli-1);
        vi=i*sr/N;
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

figure(12);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
contourf(Xm,Ym,real(Pm),25);
%caxis([0 15]);
colorbar;
%title(['window start at' num2str(tts1) 'sec']);

