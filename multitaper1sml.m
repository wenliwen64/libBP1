%close all;
clear all;


figure(1);
%M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy ];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
%r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=20;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2821044o4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+150;
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
start(6)=start(6)+50;

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
% x=x1+x2;
% 
% 
%   ts1=0.2;
%   uxx=-0.04;
%   uyy=-0.03;
%   sd=zeros(1,nel);
%   sd=-[uxx uyy ]*r';
%   
%   uxx1=-0.16;
%   uyy1=-0.15;
%   sd1=zeros(1,nel);
%   sd1=-[uxx1 uyy1 ]*r';
%   
%   y=zeros(nel,window*sr);
%  figure(14)
% for j=1:nel
%     y(j,1:window*sr)=1.2*x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1)+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr)+0.030*rand(1,window*sr);
%     subplot(nel,1,j);
%     plot(y(j,:));
% end
%coe=[    1.0000    0.9413    0.6623    0.6512    0.7380    0.7968    0.8275    0.8153    0.7940    0.5561    0.7209    0.6309];
coe291=[    1.0000    0.9241    0.6798    0.7913    0.8565    1.0021    1.3457    1.3643    1.4640    0.8563    0.8432    0.9188];%2910557b4
coe282=[    1.0000    0.8753    0.7198    0.7668    0.8069    0.8918    1.0063    0.6925    0.7285    0.5605    0.5551    0.5343];
for j=1:nel
    x(j,:)=coe282(j)*x1(j,:)+coe291(j)*x2(j,:);
end
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
figure(16);



for tl=1.5:0.1:1.5
clear X;

th=tl+2;
fl=10;
fh=20;
fran=1;
% s=linspace(0,sr,0.5*sr);


for j=1:nel
   
    subplot(nel,1,j);
    plot(tt,x(j,:),tt,coe282(j)*x1(j,:),'r',tt,coe291(j)*x2(j,:),'g');
    xlim([tl th]);
    ylim([-0.02 0.02]);
end


for i=1:nel
    X(i,:)=fft(y(i,tl*sr:th*sr-1));
    
end
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

ps=21;
qs=21;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs,ps,qs);
ux=linspace(-0.23,-0.02,ps);
uy=linspace(-0.23,-0.02,qs);


% t(-0.2172, -0.1230,'white.');
% plot(-0.06,-0.215,'white.');
% plot(-0.2188,-0.0927,'white.');



for q=1:qs
    Xm(:,q)=ux;
end
for p=1:ps
    Ym(p,:)=uy';
end


kkk=0;
%w=linspace(0,sr,window*sr);
for i=fli:fran:fhi
    kkk=kkk+1;
    clear Uv A Un a wi;
    Rxx=zeros(nel,nel);
    %Rxx=abs(X(:,i))*abs(X(:,i))';
        %Rxx=X(:,i)*X(:,i)';
   % R21a(i-fli+1)=Rxx(1,1);
%     Rxx=zeros(nel,nel);
    for j=1:nel
        for k=1:nel
        Rxx(j,k)=S(i,j,k);
        end
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
A=zeros(nel,M);
for p=1:ps
    for q=1:qs
        
            A(:,1)=v(r,kv(s(i),ux(p),uy(q)));
        
        for pp=1:ps
            for qq=1:qs
                
                A(:,2)=v(r,kv(s(i),ux(pp),uy(qq)));
                Ap=(A'*A)\A';
                PA=A*Ap;
                PAp=eye(nel)-PA;
                %sigma=1/(nel-M)*trace(PAp);
                %Psml=Ap*(Rxx-sigma*eye(nel))*Ap';
                %Pm(p,q,pp,qq)=Pm(p,q,pp,qq)+log(det(A*Psml*A'+sigma*eye(nel)));
                Pm(p,q,pp,qq)=Pm(p,q,pp,qq)+trace(PAp*Rxx);
            end
        end
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
xx=zeros(1,4);
%maxpm=max(max(max(max(Pm))));
minpm=min(min(min(min(Pm))));
for p=1:ps
    for q=1:qs
        for pp=1:ps
            for qq=1:qs
                if Pm(p,q,pp,qq)==minpm
                xx=[ux(p) uy(q) ux(pp) uy(qq)];
                end
            end
        end
    end
end


% contourf(Xm,Ym,real(Pm),15);
% %caxis([0 15]);
% colorbar;
% hold on;
% plot(-0.186,-0.157,'white.');
% plot(-0.06,-0.215,'white.');
% hold off;
% %title(['window start at' num2str(tts1) 'sec']);
end
figure(5);
hold on;
 
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);