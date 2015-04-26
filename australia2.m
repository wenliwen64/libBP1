









% clear all;
% %close all;
%  load synmag1.dat;
%  load d1d2;
%  
%  %x=d1;
%  load x1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%close all;
%figure(1);
%close all;

 %name=['78123001';'71042612';'76101000';'dpsn01wt';'dpsn01wA';'dpsn01wB';'epll01wt';'epll01wA';'epll01wB'; '72081900';'elt31   ';]
INPUT_DIR='/export/raid1/lsmeng/';



    
name=[
'MR232A';
'MR237 ';
'MR237A';
'MR237B';
'MR237C';
'MR237D';
'MR238A';
'MR238B';
'MR239A';
'MR239B';
'MR239X';
'MR239Y';
'MR240A';
'MR240B';
'MR242A';
'MR242B';
'MR243A';
'MR243B';
'MR246A';
'MR246B';
'MR247A';
'MR247B';
'MR250A';
'MR250B';
'MR253A';
'MR253B';
'MR254A';
'MR254B';
'MR256A';
'MR256B';
'MR256C';
'MR257A';
'MR257B';
'MR258A';
'MR258B';
'MR259A';
'MR259B';]

%theta=[91 132 140 140 160 160 134];
theta=[0 0 0 0 0 0 0];
alpha=10;
width=30;
%theta=[100 100 95 95 90 90 80];
%while feof(flist)==0
for ii=25:25
    figure(ii);
        clear A lonA latA val x m3 m4 m5 c2 c3 x1 x2 ddl lata lona;
%     temp1=fgets(flist);
%     [temp2 temp3]=strtok(temp1);
%     [legs temp4]=strtok(temp3);
    
  
    A=load ([INPUT_DIR strtrim(name(ii,:))]);
   
%     if length(A)<=4
%         continue;
%     end
    if(A(1,1)<A(2,1))
        latA=wrev(A(:,1));
        lonA=wrev(A(:,2));
        val=wrev(A(:,3));
    else
        latA=A(:,1);
        lonA=A(:,2);
        val=A(:,3);
    end
    disA=lonA*0;
    for kk=2:length(A)
        disA(kk)=disA(kk-1)+((lonA(kk)-lonA(kk-1))^2+(latA(kk)-latA(kk-1))^2)^0.5*110;
    end
   j=1;
    for i=1:length(A)
        %if(latA(i)>-49&&latA(i)<0&&latA(i)<-46)
      %if(latA(i)>-54.5&&latA(i)<0&&latA(i)<-52)
      if(latA(i)<0&&latA(i)>-46)
          x(j)=val(i);
          lata(j)=latA(i);
          lona(j)=lonA(i);
          j=j+1;
      end
    end
    if(j==1)
        continue;
    end
    
    % 
%     
%     
%     kk=0;
%    
%     while ~feof(fids)
%          kk=kk+1;
%         clear A ;
% %         if flag==0
% %             flag=1;
% %             B = fscanf(fids,'%d %d %d %d %d %f %f %f %f',9);
% %             C = fscanf(fids,'%s',1);
% %             continue;
% %         end
%         A = fscanf(fids,'%72c',1);
%         iold=i;
%         
%         i=str2num(A(3:5));
%         lat(kk)=str2num(A(32:39));
%         lon(kk)=str2num(A(41:49));
%         res(kk)=str2num(A(51:55));
%         
%         if(iold~=i)
%         D=[lat lon res];
%         save([DIR 'MR' num2tr(iold)],D);
%         kk=0;
%         clear lat lon res;
%         end
%         
%     end
%     fclose(fids);
    
    
    
    %%%%%%%%
%x=interp1([1:length(val)],val,[1:0.25:length(val)]);


 x=x-mean(x);
 % x=real(hilbert(x)*exp(-1i*0/180*pi));
 %aaa= synmag1;
  %x=aaa(:,5);
%     ttt=linspace(0,2*pi,500);
%     x=25*sin(10*ttt);
%     xx1=25*sin(80*ttt);
xxx=[1 2 4 10 20 40]
%x=interp1([1:length(x)],x,[1:0.25:length(x)]);
%  x=real(hilbert(x1)*exp(-1i*0/180*pi));
%  xx1=real(hilbert(xx1)*exp(-1i*50/180*pi));
%  x=x+xx1;
 %x=interp1([1:length(x)],x,[1:0.5:length(x)]);
for iii=1:1

 length1=length(x);
theta=-5;
 %x1=real(hilbert(x)*exp(1i*theta/180*pi));
%  figure(1);
%  plot(x1);
 
scale=10;
m1=-xxx(iii);
m2=xxx(iii);
k=1;
% for i=1:length1-1
%     ddl(i)=x(i+1)-x(i);
% end
for i=1+width/2:length1-width/2
    c3(k)=0;
    c2(k)=0;
    xm=mean(x(i-width/2:i+width/2));
    for j=i-width/2:i+width/2
        c3(k)=c3(k)+(j-i)^3*sign(x(j))*abs(x(j))^(3/2);
        c2(k)=c2(k)+(j-i)^2*x(j);
    end
    c3(k)=c3(k)/width;
    c2(k)=c2(k)/width;
    m3(k)=c3(k)/(abs(c2(k))^(1/2)*sign(c2(k)))^3;
    if abs(c2(k))<0.01
        c2(k)=0.01*sign(c2(k));
    end
    k=k+1;
  
end
std(m3)
%alpha=300*10*300/std(m3)^2
alpha=30*(1-atan(1.5e-4*std(m3))/(pi/2))
x2=abs(hilbert(x));
for i=1:length1-width
    if m3(i)*scale>m2||m3(i)*scale<m1
    m4(i)=0;
    else
        m4(i)=m3(i);
    end
end
%m5(1)=sum(m3);

subplot(3,1,1);
%plot((1+width/2:length-width/2),abs(m3/2),'r',1:length,x/5,'b',1:length,x2/5,'g')
%plot(1:length,x/5,'b',1:length-1,ddl,'g')
plot((1+width/2:length1-width/2),atan(abs(m3)*alpha)*scale,'r-',1:length1,x/5,'b-',1:length1,0,'g')
%ylim([m1 m2]);
for theta=0:180
 x1=real(hilbert(x)*exp(-1i*theta/180*pi));
 clear m3;
k=1;
for i=1:length1-1
    ddl(i)=x1(i+1)-x1(i);
end
for i=1+width/2:length1-width/2
    c3(k)=0;
    c2(k)=0;
    xm=mean(x1(i-width/2:i+width/2));
    for j=i-width/2:i+width/2
        c3(k)=c3(k)+(j-i)^3*sign(x1(j))*abs(x1(j))^(3/2);
        c2(k)=c2(k)+(j-i)^2*x1(j);
    end
    c3(k)=c3(k)/width;
    c2(k)=c2(k)/width;
    m3(k)=c3(k)/(abs(c2(k))^(1/2)*sign(c2(k)))^3;
    if abs(c2(k))<0.01
        c2(k)=0.01*sign(c2(k));
    end
    k=k+1;
  
end
x2=abs(hilbert(x1));
for i=1:length1-width
    if m3(i)*scale>m2||m3(i)*scale<m1
    m4(i)=0;
    else
        m4(i)=m3(i);
    end
end
m5(theta+1)=sum(atan(abs(m3)*alpha)*scale);
% figure(2);
% subplot(2,1,2);
% plot((1+width/2:length-width/2),abs(m3/2),'r',1:length,x/5,'b',1:length,x2/5,'g')
% plot(1:length,x/5,'b',1:length-1,ddl,'g')
% plot((1+width/2:length-width/2),m3*scale,'r*-',1:length,x1/5,'b',1:length,0,'g')
% ylim([m1 m2]);
end

%figure(2);
subplot(3,1,2);
plot(1:181,m5,'-')
title(name(ii,:));
subplot(3,1,3);
for theta=0:180
    if m5(theta+1)==min(m5)
        break;
    end
end
% if theta<100
%     theta=theta+180;
% end
x1=real(hilbert(x)*exp(-1i*theta/180*pi));
plot(1:length1,x1/5,'b-',1:length1,0,'g')
title(['theta=' num2str(theta)]);


%ylim([-1e5 1e5])
if iii==4
    save m5the m5;
end
end
figure(38);
scale1=4000;
hold on;
    
            plot(-lata(1:length(lata)),x1/scale1+lona(1),'b');
            %text(-lata(length(lata)),lona(1),[name(ii,:) '-' num2str(theta)]);
                        text(-lata(1),lona(1),[name(ii,:) '-' num2str(theta)]);
            xlim([ 42 56 ]);
            ylim([122 132]);






end
