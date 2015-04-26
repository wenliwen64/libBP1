%close all;
clear all;


figure(1);
M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/0.70^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=13.0;
t=linspace(0,dur,dur*sr);
x0=zeros(nel,dur*sr);
k=1;
filename='2751849m4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
start(6)=start(6)+50;
window=12.5;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);2821044o4.p
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);2751855j4.p
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x0(k,:) pr]=load_bbdata([filename '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x0(k,:) pr]=load_bbdata([filename num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x0(k,:)=x0(k,:)-mean(x0(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end
x2=zeros(nel,dur*sr);
k=1;
filename1='3120257p4.p';%'3130000h4.p';%'3120257p4.p';%'3091728j4.p';'3120257p4.p';
[xtext pr]=load_bbdata([filename1 '11' ],[],dur);

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

x1=1e4*x2;
figure(1);
plot(t,x2(1,:),t,x0(1,:),'r');
for ii=1:1
L=1.0;
ta1=2.10;
ta2=ta1+L;
xsum=zeros(1,L*sr);
yr=zeros(nel,L*sr);
xrr=zeros(nel,L*sr);
xorin=zeros(nel,L*sr);
tt=linspace(ta1,ta2,L*sr);
%sd=[   -0.0033   -0.0260   -0.0421   -0.0968   -0.1008   -0.1022   -0.1543   -0.1823   -0.2117   -0.1556   -0.1356   -0.1262];
 %sd=  [-0.0030   -0.0270   -0.0471   -0.0952   -0.1012   -0.1032   -0.1553   -0.1813   -0.2154   -0.1533   -0.1392   -0.1312]
 %sd=[     -0.0033   -0.0234   -0.0234   -0.0848   -0.0928   -0.0888   -0.1422   -0.1730   -0.2170   -0.1609   -0.1716   -0.1663]
  %sd=[     -0.0047   -0.0220   -0.0354   -0.0875   -0.0955   -0.0955   -0.1449   -0.1649   -0.1917   -0.1409   -0.1195   -0.1048] % 319
  %sd=[   -0.0047   -0.0274   -0.0421   -0.0955   -0.1022   -0.1035   -0.1556   -0.1810   -0.2130   -0.1529   -0.1359   -0.1202]%282
  
  %sd=[   -0.0047   -0.0234   -0.0247   -0.0848   -0.0968   -0.0942   -0.1449   -0.1756   -0.2184   -0.1649   -0.1716   -0.1663]%291
    % sd=[-0.0047   -0.0167   -0.0020   -0.0474   -0.0528   -0.0487   -0.0768   -0.1102   -0.1503   -0.1075   -0.1369   -0.1409];%274
     %sd=[ -0.0047   -0.0180   -0.0020   -0.0648   -0.0755   -0.0728   -0.1209   -0.1623   -0.2104   -0.1409  -0.1649   -0.1663];%2751849m4.p
     %sd=[   -0.0047   -0.0140    0.0007   -0.0274   -0.0327   -0.0274   -0.0621   -0.0875   -0.1222   -0.0888   -0.1155   -0.1222];%3091728j4
    %sd=[ -0.0033   -0.0220   -0.0073   -0.0715   -0.0795   -0.0755   -0.1235   -0.1636   -0.2096   -0.1449  -0.1636   -0.1703];%2751855j4.p
    %sd=[   -0.0033   -0.0220   -0.0087   -0.0661   -0.0755   -0.0715   -0.1169   -0.1556   -0.1957   -0.1409   -0.1569   -0.1623];%2751329c4.p
    sd=[   -0.0033   -0.0234   -0.0341   -0.0902   -0.0982   -0.0968   -0.1462   -0.1716   -0.2050   -0.1476   -0.1289   -0.1169];%3120257p4.p
    %sd=[   -0.0033   -0.0220   -0.0354   -0.0875   -0.0968   -0.0968   -0.1476   -0.1703   -0.2023   -0.1462   -0.1249   -0.1129];%3130000h4.p
  %coe1=zeros(20,nel);
  %for i=1:20
      for j=1:nel
          %
          coe1(j)=1;%min(x1(1,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+0.05*sr-1))/min(x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+0.05*sr-1));
          %end
      end
  %end
%coe=[    1.0000    0.9413    0.6623    0.6512    0.7380    0.7968    0.8275    0.8153    0.7940    0.5561    0.7209    0.6309]
    %coe=[1.0000    0.8698    0.5305    0.5989    0.6434    0.7179    1.0078    1.1782    1.1675    0.5951    0.9969    0.9291]1855
    %coe=[1.0000    0.8118    0.4532    0.4693    0.5213    0.5731    0.8778    1.0656    1.0024    0.5032    0.7813    0.7569]1849
coe=ones(1,nel);
for j=1:nel
    xsum=xsum+coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    xorin(j,:)=coe(j)*x1(j,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+L*sr-1);
end
for j=1:nel
    xr(j,:)=coe(j)*x1(j,:);
end
xsum=xsum/nel;
xres=zeros(1,L*sr);
%xr=x1;
for j=1:nel
    xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)=xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)-xsum;
end
for j=1:nel
    xres(j,:)=xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1);
end
for j=1:nel
    yr(j,:)=fft(xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1));
end
for j=1:nel
    ph=2*pi*rand(1,L*sr);
    xxr(j,:)=ifft(abs(yr(j,:)).*(cos(ph)+1i*sin(ph)));
end
for j=1:nel
    xr(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1)=xsum+xxr(j,:);
end
x=xr;

figure(13);
%plot(tt,xsum);
for j=1:nel
subplot(nel,1,j);
%ttt=linspace(1.5,2.5,sr+1);
plot(tt, coe1(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');
%plot(tt,xres(j,:),'r',tt, coe(j)*x1(j,round((ta1+sd(j))*sr):round((ta1+sd(j))*sr)+L*sr-1),'b',tt,xsum,'g')%,tt,xsum,'b',tt,x1(j,round((ta1+sd(j))*sr):round(ta1+sd(j))*sr)+L*sr-1),'g');
%xlim([1.7 3.3]);
ylimmin=min(x1(1,round((ta1+sd(1))*sr):round((ta1+sd(1))*sr)+0.2*sr-1));
%ylim([ylimmin -ylimmin]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts1=0.1;
  uxx=-0.6;
  uyy=0.3;
 % sd=-[uxx uyy (1/0.7^2-uxx^2-uyy^2)^0.5]*r';
 sd=zeros(1,nel);
for j=1:nel
    %y(j,1:window*sr)=x(j,round((ts1+sd(j))*sr):round((ts1+sd(j))*sr)+window*sr-1);%+x2(5,round((ts1+sd1(j))*sr)+1:round((ts1+sd1(j))*sr)+window*sr);
end
y=x;
     Rxx=zeros(nel,nel);
for i=1:length(tt)
    Rxx=Rxx+xorin(:,i)*xorin(:,i)';
end
for i=1:nel
    for j=1:nel
        Rxx1(i,j)=Rxx(i,j)/sqrt(sum(xorin(i,:).^2)*sum(xorin(j,:).^2));
    end
end
% figure(15);
% tt=linspace(0,dur,sr*dur)
% plot(tt,x1(1,:),'r',tt,x2(1,:),'b');
% 
fff=linspace(0,200,length(tt));
figure(11);
%xorin=xres;
zxorin=zeros(nel,length(fff));
[s c ph ci phi]=cmtm1(xorin(1,:),xorin(1,:),0.005,1.5,0,0,0);
S=zeros(nel,length(fff));
for i=1:nel
    
          [c ph ci ]=cmtm2(xorin(i,:),xorin(i,:),2);
          for k=1:length(s)
              S(i,k)=c(k);
          end
   
end

for i=1:7
    zxorin(i,:)=log10(S(i,:)./S(9,:));
   
    subplot(7,1,i);
    hold on;
    plot(fff,zxorin(i,:),'b'); 
    xlim([0 50]);
   % ylim([-3 3]);
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
% 
% 
% for tl=1.0:0.1:1.0
% clear X;
% 
% th=tl+0.5;
% fl=3;
% fh=20;
% % s=linspace(0,sr,0.5*sr);
% for i=1:nel
%     X(i,:)=fft(y(i,tl*sr:th*sr-1));
%     
% end
%  s1=linspace(0,sr,0.5*sr);
%           [s c ph ci phi]=cmtm1(y(1,tl*sr:th*sr-1),y(2,tl*sr:th*sr-1),0.005,1.5,0,0,0);
%      
%          % fspec=linspace(0,sr/2,(th-tl)*sr/2);
%           fli=round(interp1(s,1:length(s),fl));
%           fhi=round(interp1(s,1:length(s),fh));
%           S=zeros(length(s),nel,nel);
% for i=1:nel
%     for j=1:nel
%           [c ph ci ]=cmtm2(y(i,tl*sr:th*sr-1),y(j,tl*sr:th*sr-1),2);
%           for k=1:length(s)
%               S(k,i,j)=c(k);
%           end
%     end
% end
% 
% ps=200;
% qs=200;
% Xm=zeros(ps,qs);
% Ym=zeros(ps,qs);
% Pm=zeros(ps,qs);
% ux=linspace(-0.7,0.7,ps);
% uy=linspace(-0.7,0.7,qs);
% 
% for q=1:qs
%     Xm(:,q)=ux;
% end
% for p=1:ps
%     Ym(p,:)=uy';
% end
% 
% 
% 
% %w=linspace(0,sr,window*sr);
% for i=fli:fhi
%     clear Uv A Un a wi;
%     Rxx=zeros(nel,nel);
%     %Rxx=abs(X(:,i))*abs(X(:,i))';
%        % Rxx=X(:,i)*X(:,i)';
%    % R21a(i-fli+1)=Rxx(1,1);
% %     Rxx=zeros(nel,nel);
%     for j=1:nel
%         for k=1:nel
%         Rxx(j,k)=S(i,j,k);
%         end
%     end
%     %R21b(i-fli+1)=Rxx(1,1);
%     [Uv,A]=eig(Rxx);
%     
%     un=zeros(nel,nel);
%     un(:,1:nel-2)=Uv(:,1:nel-2);
%     Un=un*un';
%     %wi=fspec(i+fli-1);
%         vi=s1(i);
%         %wi=s1(i);
%         wi=1;
%     for p=1:ps
%         for q=1:qs
%             
%             a=v(r,kv(vi,ux(p),uy(q)));
%             Pm(p,q)=Pm(p,q)+wi*(a'*Un*a)/(a'*a);
%         end
%     end
% end
% % figure(13)
% % subplot(2,1,1);
% % plot(s1(fli:fhi),R21a,'b*');
% % subplot(2,1,2);
% % plot(s(fli:fhi),R21b,'b*');
% 
% figure(20);
% %contourf(X,Y,real(Pm),0.5:0.1:3.5);
% %caxis([0.5 3.5]);
% contourf(Xm,Ym,real(Pm),25);
% %caxis([0 15]);
% colorbar;
% %title(['window start at' num2str(tts1) 'sec']);
% figure(13);
% clear minpm;
% minpm=min(min(Pm));
% for p=1:ps
%     for q=1:qs
%         if Pm(p,q)==minpm
%             pm(ii)=ux(p);
%             qm(ii)=uy(q);
%             break;
%         end
%         
%     end
% end
% 
% 
% 
% 
 end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% figure(17)
% plot(pm,qm,'.');
% %xlim([-0.3 0.0]);
% %ylim([-0.3 0.0]);