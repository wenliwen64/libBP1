%close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.7^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(1i*(-r*k'));%steering vector
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
filename1='2721715j1.p';
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
x=zeros(nel,window*sr);
%x=x1;
window=5.5;

 ts1=0.1;
  uxx=0.6;
  uyy=-0.3;
  sd=-[uxx uyy (1/0.7^2-uxx^2-uyy^2)^0.5]*r';
 sd=zeros(1,nel);
for j=1:nel
    y(j,1:window*sr)=x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
x=y;
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
% %w=mywigner(x(1,:)');
% % tfrwv([x(1,:)'],t);
% % [w ttff fftt ]=tfrwv([x(1,:)',x(12,:)'],t);
% pcolor(Yfw,Xtw,abs(w));
% shading interp;
% xlim([0 20]);
figure(1);
spectrogram(x(1,:),128,120,256,200); 
Stest=spectrogram(x1(1,:),128,120,256,200); 
sizespec=size(Stest);
tspec=linspace(0,window,sizespec(2));
fspec=linspace(0,sr/2,sizespec(1));
S=zeros(sizespec(2),nel,sizespec(1));
tl=1.8;
th=2.8;
fl=5;
fh=12;
tli=round(interp1(tspec,1:sizespec(2),tl));
thi=round(interp1(tspec,1:sizespec(2),th));
fli=round(interp1(fspec,1:sizespec(1),fl));
fhi=round(interp1(fspec,1:sizespec(1),fh));
for i=1:nel
spectrogram(x(i,:),128,120,256,200); 
Sin=spectrogram(x(i,:),128,120,256,200); 
for j=tli:thi
    for k=fli:fhi
        S(j-tli+1,i,k-fli+1)=Sin(k,j);
    end
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
for i=1:fhi-fli+1
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    for j=1:thi-tli+1
    Rxx=Rxx+S(j,:,i)'*S(j,:,i);
    end
    [Uv,A]=eig(Rxx);
    
    un=zeros(nel,nel);
    un(:,1:nel-2)=Uv(:,1:nel-2);
    Un=un*un';
    wi=fspec(i+fli-1);
    for p=1:ps
        for q=1:qs
            
            a=v(r,kv(wi,ux(p),uy(q)));
            Pm(p,q)=Pm(p,q)+(a'*Un*a)/(a'*a);
        end
    end
end

figure(11);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
contourf(Xm,Ym,real(Pm),30);
%caxis([0 15]);
colorbar;
%title(['window start at' num2str(tts1) 'sec']);
