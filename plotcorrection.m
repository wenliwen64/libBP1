clear;
close all;
load ANBP0207;
load CNBP0207;
 A=ANBP0207;
 C=CNBP0207;
 
    lonA=A(:,1);
    latA=A(:,2);
     valA=A(:,3);
     disA=lonA*0;
     disA(1)=-13;
    for k=2:length(A)
        disA(k)=disA(k-1)+((lonA(k)-lonA(k-1))^2+(latA(k)-latA(k-1))^2)^0.5*110;
        if latA(k)>-0.5
            if latA(k)<0.5
                disA(k)
            end
        end
    end
     lonC=C(:,1);
    latC=C(:,2);
     valC=C(:,3);
     disC=lonC*0;
     valB=lonA*0;
    for k=2:length(C)
        disC(k)=disC(k-1)+((lonC(k)-lonC(k-1))^2+(latC(k)-latC(k-1))^2)^0.5*110;
        
    end
    for k=1:29000
         valB(k)=valC(k)-valA(k);
    end
    hold on;
    plot(latC,valC,'b',latA,valA+119,'r');
%      plot(disA,valB,'g');
    xlim([-2 6]);
    title('CM4 correction in the equatorial region');           
    ylabel('nT');
    xlabel('latitude');
    legend('observed','corrected')