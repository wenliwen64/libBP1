
clear;
close all;
DIR='~/Ge161/cf/phaseshift/';
INPUT_DIR='~/Ge161/cf/phaseshift/original/'
LIST='~/Ge161/cf/data/'
flist=fopen([LIST 'list.d'],'r');
% name=['78123001';'71042612';'76101000';'dpsn01wt';'epll01wt'; '72081900';'elt31   ';]
while feof(flist)==0
    temp1=fgets(flist);
    [temp2 temp3]=strtok(temp1);
    [legs temp4]=strtok(temp3);
    Ai=['A' legs];
    
    scale=40;
    A=load ([INPUT_DIR strtrim(Ai)]);
    
    if length(A)<=4
        continue;
    end
    
    lonA=A(:,1);
    latA=A(:,2);
    disA=lonA*0;
    for k=2:length(A)
        disA(k)=disA(k-1)+((lonA(k)-lonA(k-1))^2+(latA(k)-latA(k-1))^2)^0.5*110;
    end
    val=A(:,3);

    
    Bw=fft(val);
     p=0
%         figure;
         theta=0
%             hold on;
%             plot(disB,depth/10,'r');
%             plot(disA,val,'g');
            Y=Bw*exp(i*theta/180*pi);
            N=length(disA);
            for k=1:N/2-1
                Y(N-k+1)=Y(k+1)';
            end
            xx=ifft(Y);
            fio=fopen([DIR sprintf('%d/',theta) strtrim(legs) '.ttt'],'w');
            for kkk=1:length(lonA)
                fprintf(fio,'%f %f %f\n',lonA(kkk),latA(kkk),real(xx(kkk)));           
            end
            fclose(fio);
               
% %             plot(disA,real(xx)+(theta-p+5)*scale,'b');
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



