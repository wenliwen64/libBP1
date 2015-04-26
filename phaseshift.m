
clear;
close all;

% name=['78123001';'71042612';'76101000';'dpsn01wt';'epll01wt'; '72081900';'elt31   ';]
  name=['dpsn01wt';]
for kk=1:1
    Ai=['A' name(kk,:)];
    Bi=['B' name(kk,:)];
    scale=40;
    A=load (Ai);
    B=load (Bi);
    
    lonA=A(:,1);
    latA=A(:,2);
    disA=lonA*0;
    lonB=B(:,1);
    latB=B(:,2);
    disB=lonB*0;
    
    for k=2:length(A)
        disA(k)=disA(k-1)+((lonA(k)-lonA(k-1))^2+(latA(k)-latA(k-1))^2)^0.5*110;
    end
    for k=2:length(B)
        disB(k)=disB(k-1)+((lonB(k)-lonB(k-1))^2+(latB(k)-latB(k-1))^2)^0.5*110;
    end
    val=A(:,3);
    depth=B(:,3);
    
    Bw=fft(val);
    disA=3000-disA;
    disB=3000-disB;
         
     for p=0:90:270
         figure;
          hold on;
            plot(disB,depth/10,'r');
            plot(disA,val,'g');
         for theta=p:10:p+90
            
           
            Y=Bw*exp(1i*theta/180*pi);
            N=length(disA);
            for k=1:N/2-1
                Y(N-k+1)=Y(k+1)';
            end
            xx=ifft(Y);
            
               
             plot(disA,real(xx)+(theta-p+10)*scale,'b');
            temp=sprintf('theta=%d',theta);
            text(0,(theta-p+10)*scale,temp);
%             plot(disA,val,'g');
            title('phase shift profile');
            ylabel(name(kk,:));
            ylim([-700 4200]);
             %xlim([ 0 3000]);
             if theta==240
                 bbb=real(xx);
             end
         end
     end
%         tt=sprintf('%dto%d',p,p+90);
%         saveas(gcf,[name(kk,:) tt],'tif');
   
   
end
figure(5);
plot([1:length(disA)],disA);
    