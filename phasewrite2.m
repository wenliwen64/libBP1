
clear all;
close all;
figure(1);
DIR='~/Ge161/cf/phaseshift/';
INPUT_DIR='~/Ge161/cf/phaseshift/original/'
LIST='~/Ge161/cf/data/'
flist=fopen([LIST 'list.d'],'r');
 %name=['78123001';'71042612';'76101000';'dpsn01wt';'dpsn01wA';'dpsn01wB';'epll01wt';'epll01wA';'epll01wB'; '72081900';'elt31   ';]
 name=['78123001';'71042612';'dpsn01wA';'dpsn01wB';'epll01wA';'epll01wB'; 'elt31   ';]
flag1=[0 1 1 0 1 0 1];
theta=[91 132 140 140 160 160 134];
%theta=[0 0 0 0 0 0 0];
%theta=[100 100 95 95 90 90 80];
%while feof(flist)==0
for ii=1:7
    
        
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
   

    
    Bw=fft(val);
  
%         figure;
         
%             hold on;
%             plot(disB,depth/10,'r');
%             plot(disA,val,'g');
            Y=Bw*exp(-1i*theta(ii)/180*pi);
            N=length(disA);
            for k=1:N/2-1
                Y(N-k+1)=Y(k+1)';
            end
            xx=ifft(Y);
            %fio=fopen([DIR sprintf('%d/',theta) strtrim(name(ii,:)) '.ttt'],'w');
            fio=fopen([DIR 'theta/' strtrim(name(ii,:)) '.ttt'],'w');
            for kkk=1:length(lonA)
                fprintf(fio,'%f %f %f\n',lonA(kkk),latA(kkk),real(xx(kkk)));           
            end
            fclose(fio);
               hold on;
            plot(disA,real(xx)+ii*scale,'b');
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
fclose(flist);

x=val(230:392);
save x1 x;


