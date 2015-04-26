
close all;
clear all;
nel1=6;

lat0=20;
lon0=80
 rx0=distance(36,lon0,36,90.5)*pi/180*6371;
    ry0=distance(lat0,90.5,36,90.5)*pi/180*6371;
           %ord=[     1     2    15    13    16     7     9     8     6     3     4    10    14     5    12    11]
           
 kv=@(f0,ux,uy) 2*pi*f0*[ux uy ];%(1/0.3^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
           count=0;
figure(1)
hold on;
for i=1:nel1
%     if i==1||i==12||i==4||i==9||i==3 ||i==7
%         continue;
%     end
    count=count+1;
[Ztime,Zdata,ZSAChdr] = fget_sac(['~/matlab/kunlun/f' num2str(i)]);

    lon=ZSAChdr.station.stlo;
    lat=ZSAChdr.station.stla;
    rx=distance(lat,lon0,lat,lon)*pi/180*6371;
    ry=distance(lat0,lon,lat,lon)*pi/180*6371;
        ele=ZSAChdr.station.stel;
    r(count,1)=rx;
    r(count,2)=ry;
      %r(i,3)=ele;
    nm=ZSAChdr.station.kstnm;
    stanm(count,1:4)=nm(1:4);
    disst(count)=((rx-rx0)^2+(ry-ry0)^2)^0.5;
    sr=1/ZSAChdr.times.delta;
    
    x(count,1:length(Zdata))=Zdata;
    
    
    
    
    
    
    

% if i==1
%     count0=std(Zdata);
% else
%count=0.01*ord(i);
count0=disst(count)/2500*3;
%end
%[SeisData,HdrData,tnu,pobj]=READSAC('~/matlab/pinyon/big18/f1',1,'l')

plot(Ztime,Zdata+count0);
% if ord(i)==12
%     plot(Ztime,Zdata+count,'r');
% end
%  if count<0.125&&count>0.115
%      plot(Ztime,Zdata+count,'g');
%      ord(i)
%  end

%xlim([200 inf]);
% [B,A]=butter(4,[0.1/20 10/20]);
% Zdataf=filter(B,A,Zdata-mean(Zdata));
% figure(2)
% subplot(nel,1,i)
% plot(Ztime,Zdataf);
% xlim([100 inf]);
end
nel=count;
figure(4);
for i=1:nel/2;
    subplot(nel/2,1,i);
plot(Ztime,x(i,:)*1e-15);
%ylim([-1 1]);
end
figure(5);
for i=nel/2+1:nel;
    subplot(nel/2,1,i-nel/2);
plot(Ztime,x(i,:)*1e-15);
%ylim([-1 1]);
end


figure(3);

fff=linspace(0,sr,length(Zdata));
plot(fff,abs(fft(Zdata)));


figure(2);




% [Sefigure(1);
plot(r(:,1),r(:,2),'.');
hold on;
for i=1:nel
    text(r(i,1),r(i,2),[stanm(i,1:4) num2str(i)]);
end
r(:,1)=r(:,1)-r(1,1);
r(:,2)=r(:,2)-r(1,2);
%isData,HdrData,tnu,pobj]=READSAC('~/matlab/kunlun/',1,'b')
% t0=HdrData.B;
% t1=HdrData.E;
% dt=HdrData.DELTA;
% time=t0:dt:t1+dt;
% figure(1);
% plot(time,SeisData);



tt=Ztime;
figure(15);
%plot(tt,x1(1,:),'r',tt,1e4*x2(1,:),'b');

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
for tl=(18.5:0.1:0.1*(mm-1)+18.5)*100
    kt=kt+1
clear X;

th=tl+50;
fl=0.05;
fh=0.1;
% s=linspace(0,sr,0.5*sr);
% for i=1:nel
%     X(i,:)=fft(y(i,tl*sr:th*sr-1));
%     
% end
figure(20);
for kkk=1:nel
subplot(nel,1,kkk);
plot(x(kkk,tl*sr:th*sr-1)/std(x(kkk,tl*sr:th*sr-1)));
end
 s1=linspace(0,sr,50*sr);
          [s c ph ci phi]=cmtm1(x(1,tl*sr:th*sr-1),x(2,tl*sr:th*sr-1),0.025,1.5,0,0,0);
     
         % fspec=linspace(0,sr/2,(th-tl)*sr/2);
          fli=round(interp1(s,1:length(s),fl));
          fhi=round(interp1(s,1:length(s),fh));
          ff=round(interp1(s,1:length(s),[3          13    23    ]));
          ww=                            [1          1      1      ];
          S=zeros(length(s),nel,nel);
for i=1:nel
    for j=1:nel
          [c ph ci ]=cmtm2(x(i,tl*sr:th*sr-1)/std(x(i,tl*sr:th*sr-1)),x(j,tl*sr:th*sr-1)/std(x(j,tl*sr:th*sr-1)),2);
                    %[c ph ci ]=cmtm2(x(i,tl*sr:th*sr-1),x(j,tl*sr:th*sr-1),2);
          for k=1:length(s)
              S(k,i,j)=c(k);
          end
    end
end

ps=100;
qs=100;
zs=30;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs,zs);
ux=linspace(-0.7,0.7,ps);
uy=linspace(-0.7,0.7,qs);
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
                Rxx(j,k)=S(i,j,k);%/abs(S(i,j,k));
                end
            end
        %R21b(i-fli+1)=Rxx(1,1);
        [Uv,A]=eig(Rxx);
        As=zeros(nel,nel);
        un=zeros(nel,nel);
        us=zeros(nel,nel);
        M=11;
        %un(:,1:nel-rank(Rxx))=Uv(:,1:nel-rank(Rxx));
                un(:,1:nel-2)=Uv(:,1:nel-2);
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
                %for z=1:zs
                a=v(r,kv(vi,ux(p),uy(q)));
                Pm(p,q,1)=Pm(p,q,1)+(wi*(a'*a)/(a'*Un*a));
                %end
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
            Pmm(p,q)=Pm(p,q,1);
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




