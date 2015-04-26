






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
alpha=1;
width=20;
sc=10;
%theta=[100 100 95 95 90 90 80];
%while feof(flist)==0
for ii=32:32
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
        if(latA(i)>-52&&latA(i)<0&&latA(i)<-49)
      %if(latA(i)>-54.5&&latA(i)<0&&latA(i)<-52)
      %if(latA(i)<0&&latA(i)>-46)
          x0(j)=val(i);
          lata(j)=latA(i);
          lona(j)=lonA(i);
          j=j+1;
      end
    end
    if(j==1)
        continue;
    end


% clear all;
% %close all;
%  load synmag1.dat;
%  load d1d2;
%  
%  %x=d1;
%  load x1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
%figure(1);
%close all;

 %name=['78123001';'71042612';'76101000';'dpsn01wt';'dpsn01wA';'dpsn01wB';'epll01wt';'epll01wA';'epll01wB'; '72081900';'elt31   ';]
figure(2);
%    clear lata lona d1;
%clear x0;
 load d1d2;
 %x0=d2;
t1=linspace(0,1,length(x0));
t2=linspace(0,1,length(x0));
x=interp1(t1,x0,t2)/3;
x=x-mean(x) ;
 length1=length(x);


kk=1;

for theta=0:180
    clear x1 m3 c2 c3;
    count=0;
    x1=real(hilbert(x)*exp(-1i*theta/180*pi));
    k=1;
    ii=1;
    for i=1+width/2:length1-width/2
        
        c3(k)=0;
        c2(k)=0;
        
        for j=i-width/2:i+width/2
            c3(k)=c3(k)+(j-i)^3*sign(x1(j))*abs(x1(j))^(3/2);
            c2(k)=c2(k)+(j-i)^2*x(j);
        end
        c3(k)=c3(k)/width;
        c2(k)=c2(k)/width;
        m3(k)=c3(k)/(abs(c2(k))^(1/2)*sign(c2(k)))^3;
        if abs(c2(k))<0.01
            c2(k)=0.01*sign(c2(k));
        end
        %    if m3(k)>5
        %     m3(k)=-1;
        %   end
        
        if abs(m3(k))<=4
            m4(ii)=abs(m3(k));
            ii=ii+1;
        end
        
        if m3(k)<=0.3
            count=count+1;
        end
        k=k+1;
    end
    if mod(theta,20)==0
        
        subplot(11,2,2*kk-1);
        plot(1:length(x1),x1);
        title(['theta=' num2str(theta)])
        subplot(11,2,2*kk);
        
        hist(m4,[0:0.01:4]);
        title(['count=' num2str(count)]);
        kk=kk+1;
    end
    seq=sort(m3);
    for per=1:5
        seq=sort(m3);
        isecc(theta+1,per)=seq(round(length(m3)*per/10));
        m5(theta+1,per)=sum(atan(abs(m3)*3^per)*sc);
    end
    

 %xlim([0 4])
 %ylim([0 10])

 
end
figure(5);

 
 hold on;

 for per=1:5
 plot([0:180],(isecc(:,per)-mean(isecc(:,per)))/max(isecc(:,per)));
  plot([0:180],50*(m5(:,per)-mean(m5(:,per)))/max(m5(:,per)),'r');
 end
 hold off;
end