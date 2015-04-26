
close all;
clear all;
DIR='/export/raid1/lsmeng/';



    
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
'MR259B';
]

figure(1)
for ii=1:37
    
        
%     temp1=fgets(flist);
%     [temp2 temp3]=strtok(temp1);
%     [legs temp4]=strtok(temp3);

    
    scale=4000;
    A=load([DIR strtrim(name(ii,:))]);
   
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
    for k=2:length(A)
        disA(k)=disA(k-1)+((lonA(k)-lonA(k-1))^2+(latA(k)-latA(k-1))^2)^0.5*110;
    end
   

    hold on;
    
            plot(-latA(2:length(A)-1),val(2:length(A)-1)/scale+lonA(2),'b');
            text(-latA(length(A)-1),lonA(2),name(ii,:));
            xlim([ 42 56 ]);
            ylim([122 132]);
%             temp=sprintf('theta=%d',theta);
%             text(0,(theta-p+5)*scale,temp);
%             plot(disA,val,'g');
%             title('phase shift profile');
%             ylabel(name(kk,:));
%             ylim([-700 4200]);
%             if kk==2
%                 xlim([0 5000]);
%             end
%             if kk==4
%                 xlim([0 3000]);
%             end
%             if kk==5
%                 xlim([0 1000]);
%             end
            
%         end
        %tt=sprintf('%dto%d',p,p+90);
        %saveas(gcf,[name(kk,:) tt],'tif');
   
    %xlim([0 900]);
end






% 

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
%     fclose(fids);    MR239A     MR239Y     MR242A     MR243B     MR247A     MR250B     MR254A     MR256B     MR257B     MR259A 