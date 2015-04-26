clear all;
 gcf=figure(12);
 close(gcf);
kv=@(f0,ux,uy) 2*pi*f0*[ux uy ];% wave vector
v=@(r,k) exp(1i*(r*k'));

fc=1;
u=1;
nel=21;
lam=1/(u*fc);
r=zeros(nel,2);
r(:,1)=linspace(-0.5*lam*(nel-1)/2,0.5*lam*(nel-1)/2,nel);
% r(:,2)=linspace(-1.5,1.5,nel);
 r(:,2)=ones(1,nel)*0;


theta=40;

uxx=u*sind(theta);
uyy=u*cosd(theta);
sd=[uxx uyy ]*r';



snr=40;
sr=4;
win=128;
N=50;
dur=win*N;
t=linspace(0,dur,dur*sr);
x=zeros(nel,dur*sr);
x1=zeros(nel,dur*sr);


 y0=zeros(nel,dur*sr);
 y=zeros(nel,dur*sr);
 n=zeros(nel,dur*sr);
 x(1,:)=sin(2*pi*0.01*t+1).*sin(2*pi*fc*t);
   x1(1,:)=sin(2*pi*0.005*t+2).*sin(2*pi*fc*t);
%     x1(1,:)=sin(2*pi*fc*t);
%  figure(14);
%  hold on;
 for k=1:7
     theta1=30+k;
     uxx1=u*sind(theta1);
     uyy1=u*cosd(theta1);
     sd1=[uxx1 uyy1 ]*r';
     
     for i=1:nel
         y0(i,1:dur*sr)=specshift(x(1,:),sd(i)*sr)+specshift(x1(1,:),sd1(i)*sr);
         y(i,:)=awgn(y0(i,:),snr);
         n(i,:)=y(i,:)-y0(i,:);
         nb(i)=std(n(i,:));
         yb(i)=std(y0(i,:));
%          plot(t,y(i,:)/std(x1(1,:))+i);
     end
     
     
     
     M=2;
     kt=0;
     X=zeros(nel,win*sr);
     Rxx=zeros(nel,nel);
     for tl=0:win:(N-1)*win
         kt=kt+1;
         clear X;
         
         th=tl+win;
         %  S=zeros(tin*sr,nel,nel);
         % for ii=1:nel
         %     for j=1:nel
         %           [c ph ci ]=cmtm2(x(ii,tl*sr+1:th*sr),x(j,tl*sr+1:th*sr),3);
         %           for k=1:length(c)
         %               S(k,ii,j)=c(k);
         %           end
         %     end
         % end
         
         
         
         
         for j=1:nel
             X(j,:)=fft(y(j,tl*sr+1:th*sr)-mean(y(j,tl*sr+1:th*sr)));
         end
         Rxx=Rxx+X(:,i)*X(:,i)';
         %             for j=1:nel
         %                 for k=1:nel
         %                 Rxx(j,k)=Rxx(j,k)+S(i,j,k);
         %                 end
         %             end
         
     end
     
     ps=400;
     % qs=200;
     % Xm=zeros(ps,qs);
     % Ym=zeros(ps,qs);
     % Pm=zeros(ps,qs);
     % ux=linspace(-1,1,ps);
     % uy=linspace(-1,1,qs);
     
     % for q=1:qs
     %     Xm(:,q)=ux;
     % end
     % for p=1:ps
     %     Ym(p,:)=uy';
     % end
     phi=linspace(0,180,ps);
     Pm=zeros(1,ps);
     Pm1=zeros(1,ps);
     [Uv,A]=eig(Rxx);
     
     un=zeros(nel,nel);
     
     
     un(:,1:nel-M)=Uv(:,1:nel-M);
     
     Un=un*un';
     
     %wi=fspec(i+fli-1);
     vi=fc;
     %wi=s1(i);
     wi=1;%ww(kkk);
     for p=1:ps
         
         ux=u*sind(phi(p));
         uy=u*cosd(phi(p));
         a=v(r,kv(vi,ux,uy));
         
         Pm(p)=(wi*(a'*a)/(a'*Un*a));
         Pm1(p)=(wi*(a'*Rxx*a));%/(a'*Un*a));
         %Pm(p,q)=Pm(p,q)+(nel-2)/sum((E1(1:nel-1)-E(1:nel-1)));
         
     end
     
     figure(12);
     hold on;
     Pm=real(Pm);
     Pm1=real(Pm1);
     Pm=Pm-min(Pm);
     Pm1=Pm1-min(Pm1);
     plot(phi,Pm/max(Pm));
     plot(phi,Pm1/max(Pm1),'r');
     xlim([0 90])
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 % figure(13);
 % %contourf(X,Y,real(Pm),0.5:0.1:3.5);
 % %caxis([0.5 3.5]);
 % [ch , ch]=contourf(Xm,Ym,real(Pm),15);
 %  set(ch,'edgecolor','none');
 % %caxis([0 15]);
 % colorbar;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % ps=200;
 % qs=200;
 % Xm=zeros(ps,qs);
 % Ym=zeros(ps,qs);
 % Pm=zeros(ps,qs);
 % ux=linspace(-1,1,ps);
 % uy=linspace(-1,1,qs);
 %
 % for q=1:qs
 %     Xm(:,q)=ux;
 % end
 % for p=1:ps
 %     Ym(p,:)=uy';
 % end
 %
 %     [Uv,A]=eig(Rxx);
 %
 %     un=zeros(nel,nel);
 %
 %
 %     un(:,1:nel-M)=Uv(:,1:nel-M);
 %
 %     Un=un*un';
 %
 %     %wi=fspec(i+fli-1);
 %         vi=fc;
 %         %wi=s1(i);
 %         wi=1;%ww(kkk);
 %     for p=1:ps
 %         for q=1:qs
 %
 %             a=v(r,kv(vi,ux(p),uy(q)));
 %
 %              Pm(p,q)=(wi*(a'*a)/(a'*Un*a));
 % %              Pm(p,q)=(wi*(a'*Rxx*a));%/(a'*Un*a));
 %             %Pm(p,q)=Pm(p,q)+(nel-2)/sum((E1(1:nel-1)-E(1:nel-1)));
 %         end
 %     end
 %
 % % figure(13)
 % % subplot(2,1,1);
 % % plot(s1(fli:fhi),R21a,'b*');
 % % subplot(2,1,2);
 % % plot(s(fli:fhi),R21b,'b*');
 %
 % figure(13);
 % %contourf(X,Y,real(Pm),0.5:0.1:3.5);
 % %caxis([0.5 3.5]);
 % [ch , ch]=contourf(Xm,Ym,real(Pm),15);
 %  set(ch,'edgecolor','none');
 % %caxis([0 15]);
 % colorbar;
 
 
 

