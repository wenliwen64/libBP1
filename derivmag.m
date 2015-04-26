






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
figure(3);
%    clear lata lona d1;
%clear x0;
 load d1d2;
 %x0=d2;
t1=linspace(0,1,length(x0));
t2=linspace(0,1,length(x0));
x=interp1(t1,x0,t2)/3;
   
 length1=length(x);


k=1;

for theta=0:2:40
    clear x1;
    count=0;
 x1=real(hilbert(x)*exp(-1i*theta/180*pi));
  subplot(11,2,2*k-1);
  plot(1:length(x1),x1);
  title(['theta=' num2str(theta)])
 for i=1:length1-1
    ddl(i)=abs(x1(i+1)-x1(i));
    if ddl(i)<=3
    count=count+1;
    end
 end
 subplot(11,2,2*k);

hist(ddl,[0:.1:160]);
 title(['count=' num2str(count)]);
k=k+1;
end
end