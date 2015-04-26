%ge263Hwk2
A=1;
w=2*3.14/40;
h=5000/100;
%dt=1/100*0.707*h/4000;
dt=0.0001
p=zeros(102,102);
pp=zeros(102,102);
v=zeros(102,102);
a=zeros(102,102);
temp=zeros(102,102);
for ii=1:102
    for jj=1:102
      
       if(jj<=41||ii<=51)
            v(ii,jj)=2000;
       else
            v(ii,jj)=4000;
       end
       
       if(jj>=61)
            v(ii,jj)=4000;
       end
            a(ii,jj)=v(ii,jj)^2*dt^2/h^2;
    end
end
 figure;

for kk=1:20000
   
    for ii=1:102
        for jj=1:102
            if(ii<102&&ii>1&&jj>1&&jj<102)
            temp(ii,jj)=-p(ii,jj)+(2-4*a(ii,jj))*pp(ii,jj)+a(ii,jj)*(pp(ii-1,jj)+pp(ii+1,jj)+pp(ii,jj-1)+pp(ii,jj+1));
            end
            %free surface
            if(jj==1)
            temp(ii,jj)=0;
            end
            %bottom boundary
            if(jj==101)
            temp(ii,jj+1)=pp(ii,jj+1)+(1-a(ii,jj))/(1+a(ii,jj))*(temp(ii,jj)-pp(ii,jj+1));
            end
            %right boundary
            if(ii==101)
            temp(ii+1,jj)=pp(ii+1,jj)+(1-a(ii,jj))/(1+a(ii,jj))*(temp(ii,jj)-pp(ii+1,jj));
            end
            %right boundary
            if(ii==2)
            temp(ii-1,jj)=pp(ii-1,jj)+(1-a(ii,jj))/(1+a(ii,jj))*(temp(ii,jj)-pp(ii-1,jj));
            end
        end
    end
    if((ii>=97&&ii<=102)||(jj>=97&&jj<=102)||(ii>=1&&ii<=6)||(jj>=1&&jj<=6))
        temp(ii,jj)=0.995*temp(ii,jj);
    end
    if kk<40
    temp(50,20)= temp(50,20)+A*sin(w*kk);
    end
    p=pp;
    pp=temp;
    if mod(kk,1000)==0
        for iii=1:102
        lon(iii,:)=linspace(1,102,102);
        end
        for iii=1:102
        lat(:,iii)=linspace(1,102,102);
        end
       subplot(5,4,kk/1000)
        H=pcolor(lon, lat, temp);
        %colorbar;
        colormap(gray);
        shading interp;
    end
end

    