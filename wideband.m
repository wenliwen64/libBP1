% clear all;
% %close all;
% theta0=30*pi/180;%desired direction in rad
% fl=1.61e9;%lowest frequency in Hz
% fh=2.69e9;%highest frequency in Hz
% fl1=4.9e9;
% fh1=5.6e9;
% fc=(fl+fh)/2;
% c=3e8;%m/s
% d1=0.025;%m
% d2=d1;
% d=sqrt(d1^2*sin(theta0)^2+d2^2*cos(theta0)^2);
% rl=0.079%fl/c*d;
% rh=0.24%fh/c*d;
% rl1=0.38%fl1/c*d;
% rh1=0.49%fh1/c*d;
% % rl=fl/c*d;
% % rh=fh/c*d;
% % rl1=fl1/c*d;
% % rh1=fh1/c*d;
% N1max=40;
% N2max=40;
% N1min=10;
% N2min=10;
% % N1=(N1max-abs(theta0)/(pi/2)*(N1max-N1min));
% % N2=(N2max-abs(90-theta0/pi*180)/(90)*(N2max-N2min));
% N1=N1max;
% N2=N2max;
% phi0=atan(d1/d2*tan(theta0));
% %scale=40;
% G= @(f1,f2) (-1.8*f1^2-1.8*f2^2+0.6*sqrt(f1^2+f2^2)+0.95)*(1/16200*(atan(f1/f2))^2+1);
% Hi=zeros(N1,N2);
% f1=linspace(-0.5,0.5,N1);
% f2=linspace(-0.5,0.5,N2);
% dtheta=20*pi/180;
% alpha=d2/d1*tan(dtheta)/(-1+sqrt(1+tan(dtheta)^2+tan(dtheta)^2*tan(theta0)^2));
% Hp= @(f) sin(pi*f)/(pi*f);
% for i=1:N1
%     for j=1:N2
%         r=sqrt(f1(i)^2+f2(j)^2);
%         if (r>rl1&&r<rh1||r>rl&&r<rh)
%             H=1/G(f1(i),f2(j));
%         else
%             H=1/sqrt(10);
%         end
%         Hi(i,j)=H*Hp(alpha*(f1(i)/f2(j)-tan(phi0)));
%         R(i,j)=r;
%     end
% end
% H1=zeros(N1+1,N2+1);
% H1(2:N1+1,2:N2+1)=Hi;
%         C=real(fft2g(H1));
% % for i=1:N1max;x(i)=(i-(N1max-1)/2)*.5/((N1max-1)/2);end; %for f1-f2 axis
% % for i=1:N2max;y(i)=(i-(N2max-1)/2)*.5/((N2max-1)/2);end; %for f1-f2 axis
% figure(1)
% rotate3d on;
% %meshc(x,y,abs(QI));colormap(cool);axis([-pi pi -pi pi -.3 2]);
% surfl(abs(Hi));shading interp;colormap(pink);grid; %axis([-.5 .5 -.5 .5 -.3 2]);grid;
% %contour(x,y,abs(QI),30);colormap(cool);axis([-pi pi -pi pi -.3 1]);
% xlabel('f_2');ylabel('f_1');
% view(-20,50);
%         
% k=@(lamda,theta) 2*pi/lamda*[cos(theta) sin(theta)];% wave vector
% v=@(r,k) exp(-1i*(r*k'));%steering vector
% 
% %d=2.5;%km

% 
% r=zeros(N1*N2,2);%coordinates of the array element
% kk=1;
% for i=1:N1
%     for j=1:N2
%     r(kk,1)=(i*d1-d1)-(N1*d1-d1)/2;
%     r(kk,2)=(j*d2-d2)-(N2*d2-d2)/2;
%     W(kk)=C(i+1,j+1);
%      kk=kk+1;
%     end
% end
% 
% for fc=fc
% res=1000;
% theta3=linspace(-pi/2,pi/2,res);
% AF=zeros(1,res);
% for i=1:res
%     vv=v(r,k(c/fc,theta3(i)));
%     for j=1:N1*N2
%     AF(i)=AF(i)+W(j)*vv(j);
%     end
% end
% figure(3);
% 
% plot(theta3*180res=1000;/pi,10*log10(abs(AF).^2/(N1*N2)^2),'b');/2
% %ylim([-60 10])
% %hold on;
% end

% r=zeros(N1*N2,2);%coordinates of the array element
% kk=1;
% for i=1:N1
%     for j=1:N2
%     r(kk,1)=(i*d1-d1)-(N1*d1-d1)/2;
%     r(kk,2)=(j*d2-d2)-(N2*d2-d2)/2;
%     W(kk)=C(i+1,j+1);
%      kk=kk+1;
%     end
% end

for fc=(fl+fh)/2
res=1000;
theta3=linspace(-pi/2,pi/2,res);
AF=zeros(1,res);
for i=1:res
    %vv=v(r,k(c/fc,theta3(i)));
    for j=1:N1
        for k=1:N2
    %AF(i)=AF(i)+AIII(j,k)*exp(1i*2*pi/(c/fc)*(j*d1*cos(theta3(i))+k*d2*sin(theta3(i))));
      AF(i)=AF(i)+AIII(j,k)*exp(1i*2*pi*fc*((k-1)*cos(theta3(i))+(j-1)*sin(theta3(i))));
        end
    end
end
figure(3);

plot(theta3*180/pi,10*log10(abs(AF).^2/(N1*N2)^2),'b');
%ylim([-60 10])
%hold on;
end
%axis([-90 90 -120 -50]);
xlabel('angle');ylabel('dB');grid;
title('Directional Patterns of the Beamformer for Different Frequencies');