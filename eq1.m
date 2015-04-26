clear all;
%close all;
fid = fopen('eq1.txt','r');
kk=1;
count=0;
M1=7.0;
M2=6;
Dis=400;
Ti=60;
while ~feof(fid)
    clear A;
    A = fscanf(fid,'%405c',1);
    if size(A)==0
        break;
    end
    B(kk,9:10)=str2num(A(49:55));
    B(kk,11:16)=str2num(A(382:404));
    count=count+1;
    
    if((B(kk,9)>=M1||B(kk,10)>=M1)&&(B(kk,12)>=70||B(kk,15)>=70))
        
        
        B(kk,1)=str2num(A(6:9));
        B(kk,2)=str2num(A(11:12));
        B(kk,3)=str2num(A(14:15));
        B(kk,4)=str2num(A(17:18));
        B(kk,5)=str2num(A(20:21));
        B(kk,6)=str2num(A(23:26));
        B(kk,7)=str2num(A(28:33));
        B(kk,8)=str2num(A(35:41));
        B(kk,11:16)=str2num(A(382:404));

        B(kk,17)=count;
        kk=kk+1;
    end
end
fclose(fid);

fid = fopen('eq1.txt','r');
kk=1;
count=0;
while ~feof(fid)
    clear A;
    A = fscanf(fid,'%405c',1);
    if size(A)==0
        break;
    end
    C(kk,9:10)=str2num(A(49:55));
    count=count+1;
    
    if(C(kk,10)==0)
        C(kk,10)=(C(kk,9)-2.5)/0.63;
    end
    if(C(kk,10)>=M2)
        
        C(kk,1)=str2num(A(6:9));
        C(kk,2)=str2num(A(11:12));
        C(kk,3)=str2num(A(14:15));
        C(kk,4)=str2num(A(17:18));
        C(kk,5)=str2num(A(20:21));
        C(kk,6)=str2num(A(23:26));
        C(kk,7)=str2num(A(28:33));
        C(kk,8)=str2num(A(35:41));
        C(kk,11:16)=str2num(A(382:404));
        C(kk,17)=count;
        kk=kk+1;
    end
end
fclose(fid);
[length75 temp]=size(B);
[length6  temp]=size(C);
l=1;
 figure(3);
iii=1;
for i=1:length75
    for j=1:length6
        if C(j,17)==B(i,17)
            
            D(l,1:17)=B(i,1:17);
            D(l,18)=1;
            
            l=l+1;
            l1=l;
            D(l,:)=0;
            if j>=length6-400
                jmax=length6;
        hold on;
            else
                jmax=j+400;
            end
            for k=j+1:jmax
                
                
                Distance=110*((C(k,7)-B(i,7))^2+(C(k,8)-B(i,8))^2)^0.5;%
                Distance=distance(C(k,7),C(k,8),B(i,7),B(i,8))*110;
                Time=datenum(C(k,1),C(k,2),C(k,3),C(k,4),C(k,5),C(k,6))-datenum(B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),B(i,6));
                %if(Distance<Dis&&Time<Ti&&C(k,17)~=B(i,17)&&(D(l,9)>D(l-1,9)||D(l-1,10)>=M1))
                if(Distance<Dis&&Time<Ti&&C(k,17)~=B(i,17))  
                    l2=l;
                    D(l,1:17)=C(k,1:17);
                    D(l,18)=0;
                    D(l,19)=Distance;
                    D(l,20)=Time;
                    l=l+1;
                    D(l,:)=0;
                end
            end
            if(l1<=l2&&l2~=0)
                D(l1:l2,:)=sortrows(D(l1:l2,:),10);
                for k=l1:l2
                D(k,18)=k-l1+2;
                end
                D(l1:l2,:)=sortrows(D(l1:l2,:),20);
                plot(D(l1:l2,20),1-((1:l2-l1+1)./(l2-l1+1)),'.-','LineWidth',2);
                
                hold on;
                LD(iii)=D(l2,20);
                iii=iii+1;
            end
           
                
           
           
            l2=0;
         continue;
        end
    end
end

[lengthD  temp]=size(D);
figure(1);
subplot(3,1,1);
kk=0;
for i=1:lengthD
    
     if D(i,18)==1&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,10),'.black') ;
     end  
    if D(i,18)==2&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,10),'.r') ;
            
      
        %xlim([76 inf]);
    end
     if D(i,18)==3&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,10),'.b') ;
        
        %xlim([76 inf]);
     end
        if D(i,18)==4&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,10),'.g') ;
            
       
        %xlim([76 inf]);
        end
        if D(i,18)==5&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,10),'.yellow') ;
            
       
        %xlim([76 inf]);
        end
     if D(i,18)~=1&&D(i,1)~=0
    kk=kk+1;
    M(kk)=D(i,10); 
    end
end
 datetick('x','yy');
 ylabel('Magnitude');
subplot(3,1,2);      
%legend('red first aftershock','blue second','green third ','yellow fourth');
kk=0;
for i=1:lengthD
    
    if D(i,18)==2&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,19),'*r') ;
            
       
    end
    if D(i,18)==3&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,19),'*b') ;
            
        
    end
    if D(i,18)==4&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,19),'*g') ;
            
       
    end
    if D(i,18)==5&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,19),'*y') ;
            
      
    end
    if D(i,18)~=1&&D(i,1)~=0
    kk=kk+1;
    Dist(kk)=D(i,19); 
    end
end  
        datetick('x','yy');
        ylabel('distance from mainshock(km)');
        ylim([0 200]);
subplot(3,1,3);     
kk=0;
for i=1:lengthD
    
    if D(i,18)==2&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,20),'*r') ;
            
        
    end
    
    if D(i,18)==3&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,20),'*b') ;
            
        
    end
    
    if D(i,18)==4&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,20),'*g') ;
            
        
    end
    
    if D(i,18)==5&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,20),'*y') ;
            
       
    end
    if D(i,18)~=1&&D(i,1)~=0
    kk=kk+1;
    day(kk)=D(i,20); 
    end
end           
        
       datetick('x','yy');
        ylabel('time after mainshock(day)');
        ylim([0 10]);  
        
    figure(2);
    subplot(3,1,1);
    hist(M,100);
    %xlim([0 inf])
    xlabel('Magnitude');
    subplot(3,1,2);
    hist(Dist,200);
       xlabel('distance');
    subplot(3,1,3);
    hist(day,200);
           xlabel('time');
        
  figure(3);
  [N xx]=hist(day,200);
   for i=1:length(N)
       Num(i)=sum(N(1:i))
   end
   hold on;
   plot(xx,1-Num/sum(N),'r.-');
  figure(5);
  hold on;
  t=linspace(0,60,1000);
  Numt=zeros(1,1000);
  for i=1:1000
      for j=1:length(LD)
          if(LD(j)<=t(i))
              Numt(i)=Numt(i)+1;
          end
      end
  end
  countD=0;
  for iii=1:length(D)
      if D(iii,18)==2
      countD=countD+1;
      end
  end
  plot(t,(1-Numt/length(LD))*countD/length75,'g');
        title('chance of finding at least one large aftershock');
        xlabel('days');
       legend('Magnitude 6','Magnitude 5.5','Magnitude 5')
       ylim([0 0.8])
        
        
        