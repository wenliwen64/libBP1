close all;
clear;
%ge263Hwk2
A=1;
w=2*3.14/40;
h=5000/100; 
rou=1;
b=1/rou;
%dt=1/100*0.707*h/4000;
dt=0.0001;

txxn=zeros(102,102);
tzzn=zeros(102,102);
txzn=zeros(102,102);
txxm=zeros(102,102);
tzzm=zeros(102,102);
txzm=zeros(102,102);
vnp=zeros(102,102);
wnp=zeros(102,102);
vmp=zeros(102,102);
wmp=zeros(102,102);
vp=zeros(102,102);
vs=zeros(102,102);
lam=zeros(102,102);
miu=zeros(102,102);
aa=zeros(102,102);
ab=zeros(102,102);
temp=zeros(102,102);
for ii=1:102
    for jj=1:102
      
       if(jj<=41)
            vp(ii,jj)=2000;
       else
            vp(ii,jj)=4000;
       end
       
       
            
    end
end
 figure;
 vs=vp/3^0.5;
 miu=vs.^2*rou;
 lam=vp.^2*rou-2*miu;
 aa=vp.^2*dt^2/h^2;
 ab=vs.^2*dt^2/h^2;
for kk=1:20000
    for i=1:102
        for j=1:102
            if(i<102&&i>1&&j>1&&j<102)
                dxtxxm(i,j)=(txxm(i+1,j)-txxm(i,j))/h;
                dztxzm(i,j)=(txzm(i,j+1)-txzm(i,j))/h;
                dztzzm(i,j)=(tzzm(i,j)-tzzm(i,j-1))/h;
                dxtxzm(i,j)=(txzm(i,j)-txzm(i-1,j))/h;
    
            end
        end
    end
   
    for i=1:102
        for j=1:102
            if(i<102&&i>1&&j>1&&j<102)
                vnp(i,j)=vmp(i,j)+dt*b*(dxtxxm(i,j)+dztxzm(i,j));
                wnp(i,j)=wmp(i,j)+dt*b*(dztzzm(i,j)+dxtxzm(i,j));
                
            end
        end
    end
    
    for i=1:102
        for j=1:102
            if(i<102&&i>1&&j>1&&j<102)
                dxvx(i,j)=(vnp(i,j)-vnp(i-1,j))/h;
                dxvz(i,j)=(wnp(i+1,j)-wnp(i,j))/h;
                dzvx(i,j)=(vnp(i,j)-vnp(i,j-1))/h;
                dzvz(i,j)=(wnp(i,j+1)-wnp(i,j))/h;
            end
        end
    end
    for i=1:102
        for j=1:102
            if(i<102&&i>1&&j>1&&j<102)
                txxn(i,j)=txxm(i,j)+dt*((lam(i,j)+2*miu(i,j))*dxvx(i,j)+lam(i,j)*dzvz(i,j));
                tzzn(i,j)=tzzm(i,j)+dt*((lam(i,j)+2*miu(i,j))*dzvz(i,j)+lam(i,j)*dxvx(i,j));
                txzn(i,j)=txzm(i,j)+dt*miu(i,j)*(dzvx(i,j)+dxvz(i,j));
            end
        end
    end
   
    for i=1:102
        for j=1:102
            
                if(j==1)
                txxn(i,j)=0;
                tzzn(i,j)=0;
                txzn(i,j)=0;
                vnp(i,j)=0;
                wnp(i,j)=0;
                end
            %bottom boundary
                if(j==101)
                txxn(i,j+1)=txxm(i,j+1)+(1-aa(i,j))/(1+aa(i,j))*(txxn(i,j)-txxm(i,j+1));
                tzzn(i,j+1)=tzzm(i,j+1)+(1-aa(i,j))/(1+aa(i,j))*(tzzn(i,j)-tzzm(i,j+1));
                txzn(i,j+1)=txzm(i,j+1)+(1-ab(i,j))/(1+ab(i,j))*(txzn(i,j)-txzm(i,j+1));
                vnp(i,j+1)=vmp(i,j+1)+(1-ab(i,j))/(1+ab(i,j))*(vnp(i,j)-vmp(i,j+1));
                wnp(i,j+1)=wmp(i,j+1)+(1-aa(i,j))/(1+aa(i,j))*(wnp(i,j)-wmp(i,j+1));
                end
            %right boundary
                if(i==101)
                
                txxn(i+1,j)=txxm(i+1,j)+(1-aa(i,j))/(1+aa(i,j))*(txxn(i,j)-txxm(i+1,j));
                tzzn(i+1,j)=tzzm(i+1,j)+(1-aa(i,j))/(1+aa(i,j))*(tzzn(i,j)-tzzm(i+1,j));
                txzn(i+1,j)=txzm(i+1,j)+(1-ab(i,j))/(1+ab(i,j))*(txzn(i,j)-txzm(i+1,j));
                vnp(i+1,j)=vmp(i+1,j)+(1-aa(i,j))/(1+aa(i,j))*(vnp(i,j)-vmp(i+1,j));
                wnp(i+1,j)=wmp(i+1,j)+(1-ab(i,j))/(1+ab(i,j))*(wnp(i,j)-wmp(i+1,j));
                end
            %right boundary
                if(i==2)
                %temp(i-1,j)=pp(i-1,j)+(1-a(i,j))/(1+a(i,j))*(temp(i,j)-pp(i-1,j));
                txxn(i-1,j)=txxm(i-1,j)+(1-aa(i,j))/(1+aa(i,j))*(txxn(i,j)-txxm(i-1,j));
                tzzn(i-1,j)=tzzm(i-1,j)+(1-aa(i,j))/(1+aa(i,j))*(tzzn(i,j)-tzzm(i-1,j));
                txzn(i-1,j)=txzm(i-1,j)+(1-ab(i,j))/(1+ab(i,j))*(txzn(i,j)-txzm(i-1,j));
                vnp(i-1,j)=vmp(i-1,j)+(1-aa(i,j))/(1+aa(i,j))*(vnp(i,j)-vmp(i-1,j));
                wnp(i-1,j)=wmp(i-1,j)+(1-ab(i,j))/(1+ab(i,j))*(wnp(i,j)-wmp(i-1,j));
                end
           
        end
    end
    txxm=txxn;
    tzzm=tzzn;
    txzm=txzn;
    vmp=vnp;
    wmp=wnp;       
    if kk<40
    %txxm(50,20)=txxm(50,20)+A*sin(w*kk);
    txzm(50,20)=txzm(50,20)+A*sin(w*kk);
    
    %tzzm(50,20)=tzzm(50,20)+A*sin(w*kk);
    end  
    if mod(kk,1000)==0
        for ii=1:102
        lon(ii,:)=linspace(1,102,102);
        end
        for ii=1:102
        lat(:,ii)=linspace(1,102,102);
        end
        figure(1);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, txxm+tzzm);
        %colorbar;
        %colormap(gray);
        shading interp;
        figure(2);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, (txxm.^2+tzzm.^2).^0.5);
        %colorbar;
        %colormap(gray);
        shading interp;
        figure(3);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, vmp+wmp);
        %colorbar;
        %colormap(gray);
        shading interp;
        figure(4);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, (vmp.^2+wmp.^2).^0.5);
        %colorbar;
        %colormap(gray);
        shading interp;
        figure(5);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, txzm);
        %colorbar;
        %colormap(gray);
        shading interp;
        figure(6);
        subplot(5,4,kk/1000)
        H=pcolor(lon, lat, vmp./wmp);
        %colorbar;
        %colormap(gray);
        shading interp;
    end
 end 