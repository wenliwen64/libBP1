clear all;
close all;
load P;
% Ray tracing in a vertically heterogeneous medium with wave speed c(z)
% To shoot rays at given initial angles
% we solve the following ODE system
%   di/ds = dc/dz(z) * sin(i)
%   dx/ds = c(z)*sin(i)
%   dz/ds = c(z)*cos(i)
% where i = ray angle, counter-clockwise from Z down
%       s = travel-time
% The ODE is integrated by a midpoint method
%
% June 2008        J.P. Ampuero        ampuero@gps.caltech.edu

%I0S = [19.5:0.02:19.6];
I0S = 0.01:0.1:15;         % initial angles, in degrees
XSTA = 0; ZSTA = 0;        % initial point
%ZEND = 0;                % vertical position of final point
XEND = 120e3;                % horizontal position of final point
ZMAX = 25e3;                % stop ray tracing if the ray goes this deep
DS = 1e-2;   % timestep
resolution=100;
L=linspace(2,XEND/1e3,resolution)*1e3;
Lz=linspace(5,ZMAX/1e3,resolution)*1e3
D=zeros(1,length(L));
MM=zeros(length(I0S),length(L));
MM1=zeros(length(Lz)-1,length(L));
MM2=zeros(length(Lz),length(L));
c0=1.182e3

% example of velocity model c(z): an exponential distribution
%CFUN = @(z) 2000 + (1-exp(-z/1000))*4000;
% its derivative with respect to z
%CZFUN = @(z) exp(-z/1000)*4000/1000;
% close all
% clear all
%load P11;
P=[7.29320740145629e-60,-2.25713355056818e-54,3.15780784770224e-49,-2.64134885422995e-44,1.47242246990781e-39,-5.77410811052121e-35,1.63959228354946e-30,-3.42090336118034e-26,5.27212978329383e-22,-5.99116415754033e-18,4.98013313198646e-14,-2.98197497509444e-10,1.25141562627232e-06,-0.00349386067308002,5.89649663355069,1273.93757459125;]

%parkersfield velocity model
% 
% CFUN= @(b) P(16)+P(15)*b+P(14)*b.^2+P(13)*b.^3+P(12)*b.^4+P(11)*b.^5+(P(10)+P(9)*b+P(8)*b.^2+P(7)*b.^3+P(6)*b.^4+P(5)*b.^5+P(4)*b.^6+P(3)*b.^7+P(2)*b.^8+P(1)*b.^9).*b.^6;
% % clear all;
% % 
% % close all;
% % load b0;
% a=0:0.1:30;
% b=CFUN(a*1000);
% b(200:301)=b(200:301)+(0:101)*3.5;
% %figure(1);
% %plot(b,-a,'.');
% P=polyfit(a,b,15);
% %P=[P 0 0 0 0 0];
% load b0;
% a=0:2:30;
% %b=CFUN(a*1000);
% figure(1);
% %plot(b,-a,'.');
% P=polyfit(a,b,16);
% % a1=0:0.1:30;
% % b=polyval(P,a1);
% % figure(2);
% %plot(b,-a1);
% 
% clear b;



load b00;

a=0:0.1:30;



CFUN= @(b) P(16)+P(15)*b+P(14)*b.^2+P(13)*b.^3+P(12)*b.^4+P(11)*b.^5+(P(10)+P(9)*b+P(8)*b.^2+P(7)*b.^3+P(6)*b.^4+P(5)*b.^5+P(4)*b.^6+P(3)*b.^7+P(2)*b.^8+P(1)*b.^9).*b.^6+(b>20e3)*(b-20e3)/100*3.5;
%@(b) polyval(P,b);
%CFUN= @(b) interp1(a,b00,b)
        %aa=P(10)+P(9)*b+P(8)*b.^2+P(7)*b.^3+P(6)*b.^4+P(5)*b.^5+P(4)*b.^6+P(3)*b.^7+P(2)*b.^8+P(1)*b.^9
% its derivative with respect to z
%CZFUN = @(b) P(15)+P(14)*b*2+P(13)*b.^2*3+P(12)*b.^3*4+P(11)*b.^4*5+P(10)*b.^5*6+P(9)*b.^6*7+P(8)*b.^7*8+P(7)*b.^8*9+P(6)*b.^9*10+P(5)*b.^10*11+P(4)*b.^11*12+P(3)*b.^12*13+P(2)*b.^13*14+P(1)*b.^14*15;
CZFUN= @(b) (CFUN(b+10)-CFUN(b))/10;



% example of velocity model c(z): an exponential distribution
%CFUN = 2000 + (1-exp(-z/1000))*4000;
% its derivative with respect to zplo
%CZFUN = @(z) exp(-z/1000)*4000/1000;
%CFUN = @(z) 5000;
%CZFUN = @(z) 0;
















0.005;
% for each initial angle, shoot a ray
for k=1:length(I0S),

 % initial conditions
  inc = I0S(k)*pi/180;
  x = XSTA;
  z = ZSTA;
  clear xray zray raylength
  n=1;
  xray(n) =x;
  zray(n) =z;

 % Solve the ODE
 % Midpoint method to integrate Y'=F(Y) with step h:
 %  Y_pred = Y_n + h/2*F(Y_n)lat
 %  Y_(n+1) = Y_n + h*F(Y_pred)
 %
 % An alternative (also second order) is Heun's method:
 %  Y_pred = Y_n + h*F(Y_n)
 %  Y_(n+1) = Y_n + h*( F(Y_n) + F(Y_pred) )/2
 %
 % Or we could use ode45

  finish = 0;
  raylength(1)=0;
  while ~finish,

   % Y_n
    i_old = inc;
    x_old = x;
    z_old = z;

   % F(Y_n)
    c = CFUN(z_old);
    cz = CZFUN(z_old);
    sini = sin(i_old);
    cosi = cos(i_old);
    is = cz*sini; 
    xs = c*sini; 
    zs = c*cosi; 

   % Y_pred = Y_n + h/2*F(Y_n)
    inc = i_old +0.5*DS*is;
    x = x_old +0.5*DS*xs;
    z = z_old +0.5*DS*zs;

   % F(Y_pred)
    c = CFUN(z);
    cz = CZFUN(z);
    sini = sin(inc);
    cosi = cos(inc);
    is = cz*sini; 
    xs = c*sini; 
    zs = c*cosi; 

   % Y_(n+1) = Y_n + h*F(Y_pred)
    inc = i_old +DS*is;
    x = x_old +DS*xs;
    z = z_old +DS*zs;

    n = n+1;
    xray(n)=x;
    zray(n)=z;
    if n>=2
    raylength(n)=raylength(n-1)+((xray(n)-xray(n-1))^2+(zray(n)-zray(n-1))^2)^0.5;
    end
    %finish = z<ZEND | z>ZMAX;
    finish = x>XEND | z>ZMAX | z<ZSTA;
  end
  D=interp1(xray,zray,L,'linear',NaN);
  raylengthL=interp1(xray,raylength,L,'linear',NaN);
  MM(k,:)=D;%MM is the z coordinates for L for ray k;
  raylenthLMM(k,:)= raylengthL;
  figure(1);
  plot(xray/1e3,-zray/1e3,'-')
  hold all
end

    
for k=1:length(L)
    front=1;
    back=length(I0S);
    for jj=1:length(I0S)
        front=jj;
        if isnan(MM(jj,k))==0
            break;
        end
    end
     for jj=1:length(I0S)
         back=length(I0S)-jj+1;
        if isnan(MM(length(I0S)-jj+1,k))==0
        
            break;
        end
    end
    MM2(:,k)=interp1(MM(front:back,k),I0S(front:back),Lz,'linear','extrap')*pi/180;%MM2 take off angle for all Lz
    raylengthLLz(:,k)=interp1(MM(front:back,k),raylenthLMM(front:back,k),Lz,'linear','extrap')*pi/180;%MM2 take off angle for all Lz
end

for k=2:length(Lz)
    MM1(k-1,:)=MM2(k,:)-MM2(k-1,:);
end

hold off
axis equal
grid on


figure(3);
        iii=0;
        for xxx=0:100:ZMAX
            ccc(iii+1)=CFUN(xxx);
            xxxx(iii+1)=xxx;
            iii=iii+1;
        end
        plot(ccc/1e3,-xxxx/1e3);
        title('1d velocity model for Pakerfield');
       xlabel('(km/s) p-wave velocity');
       ylabel('(km) depth')
for iii=1:length(Lz)-1
        lat(iii,:)=L;
end

for iii=1:length(L)
        lon(:,iii)=Lz(2:length(Lz));
end
      figure(2);
      subplot(3,1,[1 2]);
      MM3=-(Lz(2)-Lz(1))./MM1*c0./cos(MM2(1:length(Lz)-1,:))/1e6.*raylengthLLz(1:length(Lz)-1,:);%slowness error
      MM13=-(Lz(2)-Lz(1))./MM1*c0./cos(MM2(1:length(Lz)-1,:))/1e6;
     % MM3=1./raylengthLLz(1:length(Lz)-1,:);%slowness error
for iii=1:length(Lz)-1


for jjj=1:length(L)
    for ijk=1:length(Lz);
        CFUNLz(ijk)=CFUN(Lz(ijk));
    end
        TOS(:,jjj)=180/pi*asin(sind(MM2(:,jjj))'.*CFUNLz/c0);% take off angle at depth Lz in degree
end
end

        [C,H]=contourf(lat/1e3, lon/1e3, MM3,20)%,1/10*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000  ]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
       
       % clabel(C,H,'fontsize',10,'color','r','rotation',0)
       title('dz/di depth uncertainty in km^2/(degree*sec)');
       xlabel('(km) distance of the array from the fault plane');
       ylabel('(km) source depth on the fault ')
        set(H,'edgecolor','none');
       set(gca,'YDir','reverse');
       %cbar =  findobj(H,'tag','Colormap');
         %colorbar;
         %caxis([1 18]);
       subplot(3,1,3);
      hold on;
      xx=linspace(0,20,5);
       ddx=ones(1,5);
       kkk=0;
       for dd=[2 5 15]
           kkk=kkk+1;
           
           plot((xx.^2+dd^2).^0.5,ddx*kkk,'^-');
           for cc=1:length(xx)
           text((xx(cc)^2+dd^2).^0.5-0.3,kkk+0.2,num2str(xx(cc)),'fontsize',8);
           end
            text((xx(cc)^2+dd^2).^0.5+2,kkk+0.2,[ 'd=' num2str(dd)],'fontsize',10);
       end
       ylim([0 5]);
       xlim([2 XEND/1e3]);
%colormap(gray);
%        figure(4);
%        xalong=linspace(0,20,resolution)*1e3;
%        MM4=zeros(length(Lz)-1,length(xalong));
%        kkk=0;
%        for dd=[ 5  5 5 ]*1e3
%            kkk=kkk+1;
%       
%         for iii=1:length(L)
%         if L(iii)>dd
%             L1(iii)=(L(iii)^2-dd^2)^0.5;
%         else
%             L1(iii)=-(dd^2-L(iii)^2)^0.5;
%         end
%         end 
%         for iii=1:length(Lz)-1
%             if kkk==1
%             MM4(iii,:)=interp1(L1,MM13(iii,:),xalong);
%             else
%             MM4(iii,:)=interp1(L1,MM3(iii,:),xalong);
%             end
%             TOSxalong(iii,:)=interp1(L1,TOS(iii,:),xalong);
%         end   
%         for iii=1:length(Lz)-1
%         lat1(iii,:)=xalong;
%         end
%         for iii=1:length(xalong)
%         lon1(:,iii)=Lz(2:length(Lz));
%         end
%         AZM=atan(dd./xalong)*180/pi;
%         for iii=1:length(Lz)-1
%             UU=RadiationPattern(0,90,0,TOSxalong(iii,:)',AZM);
%             if kkk==1
%                MM6(iii,:)=MM4(iii,:);
%             end
%             if kkk==2
%            MM6(iii,:)=(UU(:,1)');%.^0.75;
%             end
%             if kkk==3
%            MM6(iii,:)=(UU(:,2)'.^2+UU(:,3)'.^2);%.^0.375;
%             end
%            
%           % MM6(iii,:)=UU(:,1)';
%         end
%         MM5=zeros(length(Lz)-1,2*length(xalong));
%         lat2=zeros(length(Lz)-1,2*length(xalong));
%         lon2=zeros(length(Lz)-1,2*length(xalong));
%         for iii=1:length(Lz)-1
%             for jjj=1:length(xalong)
%                    lat2(iii,jjj+length(xalong))=lat1(iii,jjj);
%                    MM5(iii,jjj+length(xalong))=MM6(iii,jjj);
%                    lon2(iii,jjj+length(xalong))=lon1(iii,jjj);
%                    
%                     lat2(iii,jjj)=-lat1(iii,length(xalong)-jjj+1);
%                    MM5(iii,jjj)=MM6(iii,length(xalong)-jjj+1);
%                    lon2(iii,jjj)=lon1(iii,length(xalong)-jjj+1);
%                    
%             end
%         end
%           for iii=1:length(Lz)-1
%             for jjj=1:2*length(xalong)
%                 
%                    %if MM5(iii,jjj)>5e0
%                       % MM5(iii,jjj)=NaN;
%                    %else
%                        %MM5(iii,jjj)=log10(MM5(iii,jjj));
%                   % end
%             end
%           end
%    
%       subplot(3,1,kkk);
%          if kkk==1
%         [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,1/20*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
%        clabel(C,H,'fontsize',10,'color','r','rotation',0)
%        set(H,'edgecolor','none');
%        title(['d=' num2str(dd/1000) 'km']);
%        xlabel('(km) distance of the source along the fault plane');
%        ylabel('(km) source depth on the fault ')
%        set(gca,'YDir','reverse');
%        colorbar;
%          caxis([0 600 ]);
%          end
%             if kkk==2
%         [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,20);%*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
%        set(H,'edgecolor','none');
%         %clabel(C,H,'fontsize',10,'color','r','rotation',0)
%        title('U/R P-wave');
%        xlabel('(km) distance of the source along the fault plane');
%        ylabel('(km) source depth on the fault ')
%        set(gca,'YDir','reverse');
%        colorbar;
%         % caxis([0 800 ]);
%             end
%             if kkk==3
%         [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,20);%*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
%        set(H,'edgecolor','none');
%         %clabel(C,H,'fontsize',10,'color','r','rotation',0)
%        title('U/R S-wave');
%        xlabel('(km) distance of the source along the fault plane');
%        ylabel('(km) source depth on the fault ')
%        set(gca,'YDir','reverse');
%        colorbar;
%          %caxis([0 800 ]);
%          end
%           %caxis('manual');L1
%       %colormap(cbar);
%        end
%      
%         %shading interp;
      
        
%?????????????????????????????????????????????????
figure(5);
      ttant=(28-5)/104;
       kkk=0;
       for ddd=[ 115 120 130  ]*1e3
           kkk=kkk+1;
        xalong=linspace(0,20,resolution)*1e3;
       MM4=zeros(length(Lz)-1,length(xalong));
         
       for iii=1:length(Lz)-1
           dd=ddd-Lz(iii)/ttant;
           L1=zeros(1,length(L));
           for jjj=1:length(L)
               if L(jjj)>dd
                   L1(jjj)=(L(jjj)^2-dd^2)^0.5;
               else
                   L1(jjj)=-(dd^2-L(jjj)^2)^0.5;
               end
           end
           
          
           if kkk==1
               MM4(iii,:)=interp1(L1,MM13(iii,:),xalong);
              
               
           else
               MM4(iii,:)=interp1(L1,MM13(iii,:),xalong);
           end
           TOSxalong(iii,:)=interp1(L1,TOS(iii,:),xalong);
       end
        for iii=1:length(Lz)-1
        lat1(iii,:)=xalong;
        end
        for iii=1:length(xalong)
        lon1(:,iii)=Lz(2:length(Lz));
        end
        AZM=atan(dd./xalong)*180/pi;
        for iii=1:length(Lz)-1
            UU=RadiationPattern(0,90,0,TOSxalong(iii,:)',AZM);
            if kkk==1
               MM6(iii,:)=MM4(iii,:);
            end
            if kkk==2
           MM6(iii,:)=MM4(iii,:);%.^0.75;
            end
            if kkk==3
           MM6(iii,:)=MM4(iii,:);%.^0.375;
            end
           
          % MM6(iii,:)=UU(:,1)';
        end
        MM5=zeros(length(Lz)-1,2*length(xalong));
        lat2=zeros(length(Lz)-1,2*length(xalong));
        lon2=zeros(length(Lz)-1,2*length(xalong));
        for iii=1:length(Lz)-1
            for jjj=1:length(xalong)
                   lat2(iii,jjj+length(xalong))=lat1(iii,jjj);
                   MM5(iii,jjj+length(xalong))=MM6(iii,jjj);
                   lon2(iii,jjj+length(xalong))=lon1(iii,jjj);
                   
                    lat2(iii,jjj)=-lat1(iii,length(xalong)-jjj+1);
                   MM5(iii,jjj)=MM6(iii,length(xalong)-jjj+1);
                   lon2(iii,jjj)=lon1(iii,length(xalong)-jjj+1);
                   
            end
        end
          for iii=1:length(Lz)-1
            for jjj=1:2*length(xalong)
                
                   %if MM5(iii,jjj)>5e0
                      % MM5(iii,jjj)=NaN;
                   %else
                       %MM5(iii,jjj)=log10(MM5(iii,jjj));
                  % end
            end
          end
   
      subplot(3,1,kkk);
         if kkk==1
        [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,[0:25: 900])%1/20*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
       clabel(C,H,'fontsize',10,'color','r','rotation',0)
       set(H,'edgecolor','none');
       title(['d=' num2str(ddd/1000) 'km']);
       xlabel('(km) distance of the source along the fault plane');
       ylabel('(km) source depth on the fault ')
       set(gca,'YDir','reverse');
       colorbar;
         caxis([300 900 ]);
         end
            if kkk==2
        [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,[0:25: 900]);%*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
       set(H,'edgecolor','none');
        clabel(C,H,'fontsize',10,'color','r','rotation',0)
      title(['d=' num2str(ddd/1000) 'km']);
       xlabel('(km) distance of the source along the fault plane');
       ylabel('(km) source depth on the fault ')
       set(gca,'YDir','reverse');
       colorbar;
         caxis([300 900 ]);
            end
            if kkk==3
        [C,H]=contourf(lat2/1e3, lon2/1e3, MM5,[0:25: 900]);%*[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000  7000  8000 9000 10000 11000  12000  13000 14000  16000  18000 20000 22000]);%[50 100 300 400 500 600 700 800 900 1000 1200 1500 2000 2500 3000 ]
       set(H,'edgecolor','none');
        clabel(C,H,'fontsize',10,'color','r','rotation',0)
      title(['d=' num2str(ddd/1000) 'km']);
       xlabel('(km) distance of the source along the fault plane');
       ylabel('(km) source depth on the fault ')
       set(gca,'YDir','reverse');
       colorbar;
         caxis([300 900 ]);
         end
          %caxis('manual');
      %colormap(cbar);
          end
     
     
        %shading interp;
      
        
