%close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy,uz) 2*pi*f0*[ux uy uz];%(1/0.3^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=12;
t=linspace(0,dur,dur*sr);
x1=zeros(nel,dur*sr);
k=1;
filename='2721715j2.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+100;
window=7.5;
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
filename1='2721715j3.p';
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
  %sd=[   -0.0047   -0.0274   -0.0421   -0.0955   -0.1022   -0.1035   -0.1556   -0.1810   -0.2130   -0.1529   -0.1359   -0.1202];%282
  sd=[ 0.0068   -0.0142   -0.0324   -0.0767   -0.0804   -0.0845   -0.1164   -0.1342   -0.1581   -0.1132   -0.0850   -0.0652];%272
  %sd=[    0.0068   -0.0142   -0.0324   -0.0767   -0.0804   -0.0845   -0.1164   -0.1342   -0.1581   -0.1132   -0.0850   -0.0652]
x=zeros(nel,dur*sr);
x=x2*cosd(50)-x1*cosd(40);%+1.5*x2;
% for j=1:nel
%    % x(j,:)=coe282(j)*x2(j,:)+coe291(j)*x1(j,:);
%         x(j,:)=2*coe282(j)*x2(j,:)+coe291(j)*x1(j,:);
% end

  ts1=1;
    uxx=-0.2;
    uyy=-0.2;
sd=-[uxx uyy (1/0.3^2-uxx^2-uyy^2)^0.5 ]*r';
%   
%   
%   
%   ts2=0.4;
%   uxx1=0.2;
%   uyy1=0.3;
%   sd1=-[uxx1 uyy1 ]*r';
%   
%   
%  %ts2=0.4;
%   uxx2=0.6;
%   uyy2=0.3;
   %sd=-[uxx uyy ]*r';
%sd=zeros(1,nel);
for j=1:nel
    %y(j,1:window*sr)=x1(1,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1)+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
    y(j,1:window*sr)=x(1,round((ts1+sd(j))*sr)+1:round((ts1+sd(j))*sr)+window*sr);
end
%y=x;

figure(15);
tt=linspace(0,dur,sr*dur);
plot(tt,x1(1,:),'r',tt,1e4*x2(1,:),'b');

figure(14);
for i=1:nel
    subplot(nel,1,i)
    plot(tt,x(i,:));
end

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

mm=29;
DOA=zeros(mm,4,2);
kt=0;
for tl=3.250:0.1:0.1*(mm-1)+3.25
    kt=kt+1
clear X;

th=tl+0.5;
fl=4;
fh=15;
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
          [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1),y(j,tl*sr:th*sr-1),3);
          for k=1:length(s)
              S(k,i,j)=c(k);
          end
    end
end

ps=30;
qs=30;
zs=30;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs,zs);
ux=linspace(-0.5,0.5,ps);
uy=linspace(-0.5,0.5,qs);
uz=linspace(-0.2,2,zs);
for q=1:qs
    Xm(:,q)=ux;%-0.2403;
end
for p=1:ps
    Ym(p,:)=uy';%-0.0374;
end


kkk=0;
%w=linspace(0,sr,window*sr);
for ii=1:1
    Pm=zeros(ps,qs,zs);
    for i=fli:1:fhi
        kkk=kkk+1;
        clear Uv A Un a wi;
        Rxx=zeros(nel,nel);
        %Rxx=abs(X(:,i))*abs(X(:,i))';
        %Rxx=X(:,i)*X(:,i)';
        % R21a(i-fli+1)=Rxx(1,1);
        %     Rxx=zeros(nel,nel);
            for j=1:nel
                for k=1:nel
                Rxx(j,k)=S(i,j,k)/abs(S(i,j,k));
                end
            end
        %R21b(i-fli+1)=Rxx(1,1);
        [Uv,A]=eig(Rxx);
        As=zeros(nel,nel);
        un=zeros(nel,nel);
        us=zeros(nel,nel);
        M=11;
        %un(:,1:nel-rank(Rxx))=Uv(:,1:nel-rank(Rxx));
                un(:,1:nel-1)=Uv(:,1:nel-1);
        us(:,nel-rank(Rxx)+1:nel)=Uv(:,nel-rank(Rxx)+1:nel);
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
                for z=1:zs
                a=v(r,kv(vi,ux(p),uy(q),uz(z)));
                Pm(p,q,z)=Pm(p,q,z)+(wi*(a'*a)/(a'*Un*a));
                end
            end
        end
    end

maxpm=max(max(max(Pm)));

    
for p=1:ps
        for q=1:qs
            for z=1:zs
                if Pm(p,q,z)==maxpm
                    DOA(kt,ii,1)=ux(p);
                    DOA(kt,ii,2)=uy(q);
                    uuz=z;
                end
            end
        end
end
% if maxpm<=10;
%     DOA(kt,ii,1)=DOA(kt-1,ii,1);
%     DOA(kt,ii,2)=DOA(kt-1,ii,2);
% end
end
    
for p=1:ps
        for q=1:qs
            Pmm(p,q)=Pm(p,q,uuz);
        end
end

% figure(13)
% subplot(2,1,1);
% plot(s1(fli:fhi),R21a,'b*');
% subplot(2,1,2);
% plot(s(fli:fhi),R21b,'b*');
% 
rr=interp1([1 mm/2 mm],[1 0 0],1:mm)
gg=interp1([1 mm/2 mm],[0 1 0],1:mm)
bb=interp1([1 mm/2 mm],[0 0 1],1:mm)

% for i=1:20
%     hold on;
%     plot(i,i,'color',[r(i) g(i) b(i)]);
% end
figure(kt+20);
%contourf(X,Y,real(Pm),0.5:0.1:3.5);
%caxis([0.5 3.5]);
[ch ch]=contourf(Xm,Ym,real(Pmm),15);
 set(ch,'edgecolor','none');
 title(['uz=' num2str(uz(uuz))]);
figure(19);
hold on;
topc=14/15*(max(max(Pmm))-min(min(Pmm))+min(min(Pmm)));



% if kt<=15
 %contour(Xm,Ym,real(Pm),[topc topc],'color',[rr(kt) gg(kt) bb(kt)]);
 if kt<=15
 contour(Xm,Ym,real(Pmm),[topc topc],'r'); 
end
if kt>=16&&kt<=18
 contour(Xm,Ym,real(Pmm),[topc topc],'b'); 
end
 if kt>=19&&kt<=21
      contour(Xm,Ym,real(Pmm),[topc topc],'g');
 end
 if kt>=22&&kt<=29
           contour(Xm,Ym,real(Pmm),[topc topc],'magenta');
 end
%  end
%caxis([0 15]);
colorbar;
% hold on;
% plot(-0.186,-0.157,'white.');
% plot(-0.06,-0.215,'white.');
% hold off;
% title([num2str(tl) num2str(kt)]);
% title(['window start at' num2str(tts1) 'sec']);
end



figure(5);
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);
hold on;

plot(DOA(:,1,1)-0.4473,DOA(:,1,2) -0.0774,'b-O');
for iii=1:mm
    text(DOA(iii,1,1)-0.4473+iii*1e-4,DOA(iii,1,2)-0.0374,num2str(iii));

end

% plot(DOA(:,2,1),DOA(:,2,2),'g-O');
% plot(DOA(1:4,3,1),DOA(1:4,3,2),'b-O');
% plot(DOA(:,4,1),DOA(:,4,2),'black-O');