clear all;
fid = fopen('eq1.txt','r');
kk=1;
count=0;
M1=7.5;
M2=5.5;
Dis=1000;
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
    if(C(kk,9)>=M2||C(kk,10)>=M2)
        
        
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
for i=1:length75
    for j=1:length6
        if C(j,17)==B(i,17)
            
            D(l,1:17)=B(i,1:17);
            D(l,18)=1;
            l=l+1;
            D(l,:)=0;
            if j>=length6-50
                jmax=length6;
            else
                jmax=j+50;
            end
            for k=j+1:jmax
                %Distance=110*((C(k,7)-B(i,7))^2+(C(k,8)-B(i,8))^2)^0.5;%
                Distance=distance(C(k,7),C(k,8),B(i,7),B(i,8))*110;
                Time=datenum(C(k,1),C(k,2),C(k,3),C(k,4),C(k,5),C(k,6))-datenum(B(i,1),B(i,2),B(i,3),B(i,4),B(i,5),B(i,6));
                if(Distance<Dis&&Time<Ti&&C(k,17)~=B(i,17)&&(D(l,9)>D(l-1,9)||D(l-1,10)>=M1))
                    D(l,1:17)=C(k,1:17);
                    D(l,18)=0;
                    D(l,19)=Distance;
                    D(l,20)=Time;
                    l=l+1;
                    D(l,:)=0;
                end
            end
            
         continue;
        end
    end
end
[lengthD  temp]=size(D);
figure(1);
subplot(3,1,1);
for i=1:lengthD
    
    if D(i,18)==0&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,9),'*') ;
            
        datetick('x','yy');
        ylabel('Magnitude');
        %xlim([76 inf]);
    end
end
subplot(3,1,2);      
for i=1:lengthD
    
    if D(i,18)==0&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,19),'*') ;
            
        datetick('x','yy');
        ylabel('distance from mainshock(km)');
        ylim([0 200]);
    end
end        
subplot(3,1,3);      
for i=1:lengthD
    
    if D(i,18)==0&&D(i,1)~=0
        hold on;
       plot(datenum(D(i,1),D(i,2),D(i,3),D(i,4),D(i,5),D(i,6)),D(i,20),'*') ;
            
        datetick('x','yy');
        ylabel('time after mainshock(day)');
        ylim([0 10]);
    end
end           
        
        
        
        
        
        
        
        
        
        
        
        