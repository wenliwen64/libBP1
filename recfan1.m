%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Rectangular Array
%%
%%
close all
clear all
QIKRMAX=0;QIKRMIN=100000000000;
N1MAX=51;N1MIN=11;N2MAX=51;N2MIN=11;BW=10;
%fl=.10;fh=.48;
T0 = 30;%input('Enter main beam direction: ');
%T0=85;
% fl=1.61e9/1.2e10;%lowest frequency in Hz
% fh=2.69e9/1.2e10;%highest frequency in Hz
% fl1=4.9e9/1.2e10;
% fh1=5.6e9/1.2e10;

c=5e3;%m/s
d1=500;%m
d2=d1;
d=sqrt(d1^2*sind(T0)^2+d2^2*cosd(T0)^2);
fl=0.5/(c/d);
fh=5.0/(c/d);
fl1=2/(c/d);
fh1=5/(c/d);
% rl=0.079%fl/c*d;
% rh=0.24%fh/c*d;
% rl1=0.38%fl1/c*d;
% rh1=0.49%fh1/c*d;

rl=fl;
rh=fh;
rl1=fl1;
rh1=fh1;

N1=round(N1MAX+(abs(T0)/90)*(N1MIN-N1MAX));
N2=round(N2MAX+((90-abs(T0))/90)*(N2MIN-N2MAX));
if round(N1/2)*2==N1;N1=N1+1;end,
if round(N2/2)*2==N2;N2=N2+1;end,
NN2=N2MAX;NN1=N1MAX;N2S=NN2-1;N1S=NN1-1;
TT=tan(T0*pi/180);TB=tan(BW*pi/180);
PX=TB/(-1+sqrt(1+TB*TB+(TB*TT)^2));
if T0<0;SG=-1;else SG=1;end,

IF2=0;
for F2=-.5+1/(2*NN2):1/NN2:.5-1/(2*NN2)+.0001;
    IF2=IF2+1;IF1=N1S/2-SG+1;
    for F1=0:SG/NN1:SG*(.5-1/(2*NN1))+.0001;
        IF1=IF1+SG;
        r=sqrt(F1^2+F2^2);
        if abs(F2)<.0001;
            QI(IF1,IF2)=QI(IF1,IF2-1);
        else
            if r>rl&&r<rh
                QI(IF1,IF2)=sin(PX*pi*((F1/F2)-TT))/(PX*pi*((F1/F2)-TT))/((-1.8*F1^2-1.8*F2^2+0.6*sqrt(F1^2+F2^2)+0.95)*(1/16200*(atan(F1/F2))^2+1));
            else
                QI(IF1,IF2)=0.1/10^0.5;
            end
        end,
    end,
end,

IF2=0;
for F2=-.5+1/(2*NN2):1/NN2:.5-1/(2*NN2)+.0001;
    IF2=IF2+1;IF1=N1S/2+SG+1;
    for F1=0:-SG/NN1:-SG*(.5-1/(2*NN1));
        IF1=IF1-SG;
        QI(IF1,IF2)=QI(N1S-IF1+2,N2S-IF2+2);
    end,
end,

%for i=1:NN1;x(i)=(i-(NN1-1)/2)*pi/((NN1-1)/2);end; %for u-v axis
%for i=1:NN2;y(i)=(i-(NN2-1)/2)*pi/((NN2-1)/2);end; %for u-v axis

for i=1:NN1;x(i)=(i-(NN1-1)/2)*.5/((NN1-1)/2);end; %for f1-f2 axis
for i=1:NN2;y(i)=(i-(NN2-1)/2)*.5/((NN2-1)/2);end; %for f1-f2 axis

%for i=1:NN1;x(i)=(i-(NN1-1)/2)*pi/((NN1-1)/2);end; %for u-v axis

figure(1)
rotate3d on;
%meshc(x,y,abs(QI));colormap(cool);axis([-pi pi -pi pi -.3 2]);
surfl(x,y,abs(QI));shading interp;colormap(pink);grid; axis([-.5 .5 -.5 .5 -.3 2]);grid;
%contour(x,y,abs(QI),30);colormap(cool);axis([-pi pi -pi pi -.3 1]);
xlabel('f_2');ylabel('f_1');
view(-20,50);

mi=0;
for IN2=-N2S/2:N2S/2;
    for IN1=-N1S/2:N1S/2;
        QIKR=0;QIKI=0;
        for IF2=-N2S/2:N2S/2;
            for IF1=-N1S/2:N1S/2;
                QIKR=QIKR+QI(IF1+N1S/2+1,IF2+N2S/2+1)*cos(2*pi*((IN2*IF2)/NN2+(IN1*IF1)/NN1));
            end,
        end,
        if abs(QIKR)<1;QIKR=0;end,
        %if abs(QIKR)>QIKRMAX;QIKRMAX=abs(QIKR)
        %if abs(QIKR)<QIKRMIN;QIKRMIN=abs(QIKR)
        %end,
        %end,
        AI(IN1+N1S/2+1,IN2+N2S/2+1)=QIKR/(NN2*NN1);
        if abs(IN2)<=(N2-1)/2;if abs(IN1)<=(N1-1)/2;mi=mi+1;aii(mi)=AI(IN1+N1S/2+1,IN2+N2S/2+1);end,end,
    end,
end,
for jj=1:N2;for ii=1:N1;AIII(N1-ii+1,N2-jj+1)=aii(ii+(jj-1)*N1);end,end,

figure(2);
for ff=1:20,
    f=((fh1-fl1)/19)*ff+(20*fl1-fh1)/19;
    for t=-90:90;
        tet=t*pi/180;
        tt=t+91;
        for m=1:N1,
            for n=1:N2,
                Ht(m,n)=AIII(m,N2-n+1)*exp(j*2*pi*f*((m-1)*sin(tet)-(n-1)*cos(tet)));
            end;
        end;
        HP(tt)=sum(sum(Ht));
    end;
    AA = [-90:1:90];
    plot(AA,20*log10(abs(HP)),'b');hold on;
end;
axis([-90 90 -60 10]);
xlabel('angle');ylabel('dB');grid;
title('Directional Patterns of the Beamformer for Different Frequencies');
save aiii AIII