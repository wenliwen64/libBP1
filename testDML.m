% test the DML algorthm in the haiti case

close all;
clear all;


figure(1);
%M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy ];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,2);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
% r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;
sr = 200; %sampling rate /sec
dur=60;
t=linspace(0,dur,dur*sr);
x=zeros(nel,dur*sr);
x1=zeros(nel,dur*sr);

  uxx=-0.04;
  uyy=-0.04;
  sd=zeros(1,nel);
  sd=-[uxx uyy ]*r';
  
  uxx1=-0.16;
  uyy1=-0.16;
  sd1=zeros(1,nel);
  sd1=-[uxx1 uyy1 ]*r';

 
 y=zeros(nel,dur*sr);
 x(1,:)=sin(2*pi*0.03*t+1).*sin(2*pi*0.3*t);
  x1(1,:)=sin(2*pi*0.02*t+3).*sin(2*pi*0.3*t);
 figure(14)
 hold on;
for j=1:nel
    y(j,1:dur*sr)=specshift(x(1,:),sd(j)*sr)+specshift(x1(1,:),sd1(j)*sr);
    
    plot(y(j,:)/std(x(1,:))+j);
end


figure(15);
tt=linspace(0,dur,sr*dur);
plot(t,x(1,:),'r',t,x1(1,:),'b');



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

Nw=3;
for tl=1.5:0.1:1.5
clear X;


fl=0.3;
fh=0.3;
% s=linspace(0,sr,0.5*sr);
% for i=1:nel
%     X(i,:)=fft(y(i,tl*sr:th*sr-1));
%     
% end
 s=linspace(0,sr,length(y));
          
%           [c ph ci ]=cmtm2(xi/std(xi),xj/std(xj),Nw);
         fli=round(interp1(s,1:length(s),fl));
          fhi=round(interp1(s,1:length(s),fh));
          S=zeros(length(s),nel,nel);
for i=1:nel
    for j=1:nel
          [c ph ci ]=cmtm2(y(i,:),y(j,:),Nw);
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