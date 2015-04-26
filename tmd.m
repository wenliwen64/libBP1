






clear all;
close all;
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
for ii=32:37
   % figure(ii);
        clear A lonA latA val x m3 m4 m5 c2 c3 x1 x2 ddl lata lona;
%     temp1=fgets(flist);40
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
    clear lata lona x0;
    for i=1:length(A)
        if(latA(i)>-53.5&&latA(i)<0&&latA(i)<-47)
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
%  load x1;40

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
%figure(1);
%close all;

alpha=1;
 %name=['78123001';'71042612';'76101000';'dpsn01wt';'dpsn01wA';'dpsn01wB';'epll01wt';'epll01wA';'epll01wB'; '72081900';'elt31   ';]
figure(3);
   
%clear x0;
 load d1d2;
 %x0=d2;
t1=linspace(0,1,length(x0));
t2=linspace(0,1,length(x0)/4);
x=interp1(t1,x0,t2);
lata1=interp1(t1,lata,t2); 
lona1=interp1(t1,lona,t2);
 length1=length(x);


k=1;

for theta=0:180
    clear x1;
    count=0;
 x1=real(hilbert(x)*exp(-1i*theta/180*pi));
 
 for i=1:length1-1
    ddl(i)=abs(x1(i+1)-x1(i));
    if ddl(i)<=3
    count=count+1;
    end
 end
 if(mod(theta,20)==0)
     subplot(11,2,2*k-1);
     plot(1:length(x1),x1);
     title(['theta=' num2str(theta)])
     subplot(11,2,2*k);
     hist(ddl,[0:.1:160]);
     title(['count=' num2str(count)]);
     k=k+1;
 end
 
 for per=1:80
     seq=sort(ddl);
     isec(theta+1,per)=seq(round(length(ddl)*(0.1+per/100)));
     m5(theta+1,per)=sum(atan(abs(ddl)*3^per)*sc);
 end

 

 
 
 

end

 figure(ii);

 tot=zeros(181,1);
 hold on;

 for per=1:80
 plot([0:180],(isec(:,per)-mean(isec(:,per)))/max(isec(:,per)));
  plot([0:180],5*(m5(:,per)-mean(m5(:,per)))/max(m5(:,per)),'r');
  tot=tot+(isec(:,per)-mean(isec(:,per)))/max(isec(:,per));
 end
 
 for theta=0:180
    if tot(theta+1)==min(tot)
        themin=theta;
        break;
    end
end
for theta=0:180
    if tot(theta+1)==max(tot)
        themax=theta;
        break;
    end
end
 plot([0:180],tot/50,'g','LineWidth',4);
 hold off;
 r=tot(themax+1)-tot(themin+1);
 ki=1;
 for theta=0:180
    if abs(tot(theta+1)-(tot(themin+1)))<=0.05*r
        thek(ki)=theta;
        ki=ki+1;
    end
 end
  title(['theta=' num2str(themin) '(' num2str(min(thek)) '-' num2str(max(thek)) ')' ]);
 figure(38);
scale1=4000;
hold on;
             x1=real(hilbert(x)*exp(-1i*themin/180*pi));
            plot(-lata1(1:length(lata1)),x1/scale1+lona1(1),'b');
            %text(-lata(length(lata)),lona(1),[name(ii,:) '-' num2str(theta)]);
                        text(-lata1(1),lona(1),[name(ii,:) '-' num2str(themin) '(' num2str(min(thek)) '-' num2str(max(thek)) ')'  ]);
            xlim([ 42 56 ]);
            ylim([122 132]);

clear thek;


 
 
 
 
 
 
 
end
