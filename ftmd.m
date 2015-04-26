






clear all;
close all;
figure(1);
DIR='~/Ge161/cf/phaseshift/';
INPUT_DIR='~/Ge161/cf/phaseshift/original/'
LIST='~/Ge161/cf/data/'
flist=fopen([LIST 'list.d'],'r');
 %name=['78123001';'71042612';'76101000';'dpsn01wt';'dpsn01wA';'dpsn01wB';'epll01wt';'epll01wA';'epll01wB'; '72081900';'elt31   ';]
 %name=['78123001';'71042612';'dpsn01wA';'dpsn01wB';'epll01wA';'epll01wB'; 'elt31   ';]
 %name=['v2113_1 ';'kk831116';'8700161A';'NBP0207 ';'dsdp08gc';'v2403   '; 'erdc11wt';]
 
%  name=['v2113_1 ';'8700161A';'kk831117';'kk831118';'kk831119';'kk831120';'NBP0207 ']%;'dsdp08gc';'v2403   '; 'erdc11wt';]
% flag1=[0 0 0 1 0 1 1];

name=['dsdp08g1';'v2404   ';'erdc11w1';'aria01w1';'NBP0209 ']
flag1=[0 1 0 1 1 1];
theta=[91 132 140 140 160 160 134];
 %cutf=[1/3 0 1/3 0.2 0.1 1/3 0.6];
 %cute=[1 3/8 2/3 1 3/5 1 1];
%theta=[0 0 0 0 0 0 0];
%theta=[100 100 95 95 90 90 80];
%while feof(flist)==0
for ii=1:7
            clear A lonA latA val x m3 m4 m5 c2 c3 x1 x2 ddl lata lona disA disA1 x0;
        
%     temp1=fgets(flist);
%     [temp2 temp3]=strtok(temp1);
%     [legs temp4]=strtok(temp3);
    Ai=['A' name(ii,:)];
    
    scale=400;
    A=load ([INPUT_DIR strtrim(Ai)]);
    
    if length(A)<=4
        continue;
    end
    if(flag1(ii)~=0)
        lonA=wrev(A(:,1));
        latA=wrev(A(:,2));
        val=wrev(A(:,3));
    else
        lonA=A(:,1);
        latA=A(:,2);
        val=A(:,3);
    end
    disA=lonA*0;
    for k=2:length(A)
        disA(k)=disA(k-1)+((lonA(k)-lonA(k-1))^2+(latA(k)-latA(k-1))^2)^0.5*110;
    end

   % figure(ii);

%     temp1=fgets(flist);40
%     [temp2 temp3]=strtok(temp1);
%     [legs temp4]=strtok(temp3);
    
  
    x0=val;


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
t2=linspace(0,1,length(x0)/2);
%t2=linspace(cutf(ii),cute(ii),0.5*length(x0)*(cute(ii)-cutf(ii)));
x=interp1(t1,x0,t2);
x=x-mean(x);
disA1=interp1(t1,disA,t2);
length1=length(x);
x=real(hilbert(x)*exp(-1i*0/180*pi));
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
 
 for per=1:60
     seq=sort(ddl);
     isec(theta+1,per)=seq(round(length(ddl)*(0.2+per/100)));

 end

 

 
 
 

end

 figure(1);
 subplot(7,2,2*ii);
 tot=zeros(181,1);
 hold on;

 for per=1:60
 plot([0:180],(isec(:,per)-mean(isec(:,per)))/max(isec(:,per)));

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
 plot([0:180],tot/40,'r','LineWidth',3);
       %ylim([-0.2 0.2]);
 hold off;
 r=tot(themax+1)-tot(themin+1);
 ki=1;
 for theta=0:180
    if abs(tot(theta+1)-(tot(themin+1)))<=0.05*r
        thek(ki)=theta;
        ki=ki+1;
    end
 end

 
  title([ 'theta=' num2str(themin) '(' num2str(min(thek)) '-' num2str(max(thek)) ')' ]);
   subplot(7,2,2*ii-1);
                x1=real(hilbert(x)*exp(-1i*0/180*pi));
                plot(disA1,x1);
                if(ii==4)
                   % xlim([650 inf])
                end
                title(name(ii,:));
          
 figure(38);
scale1=400;
hold on;

            plot(disA1,x1+scale1*ii,'b');

                        text(-10,scale1*ii,[name(ii,:) '-' num2str(themin) '(' num2str(min(thek)) '-' num2str(max(thek)) ')'  ]);
           

clear thek;
end
