    function Pm=musicdeep(x0,tl,th,s,fli,fs,fhi,Nw,sr,Mm,ps,qs,zs,tlib)

    [nel Lh]=size(x0);
    y=zeros(nel,Lh);
    
    x00=x0;
    
    S=cmtmall(x00(:,tl*sr:th*sr-1),Nw);
    
    
    if strcmp(Mm,'rank')
        M=rank(Rxx);
        ['M=' num2str(M)];
    else
        M=Mm;
    end
    
        Pm=zeros(ps,zs);
        Pw=zeros(ps,zs);
    for i=fli:fs:fhi
       s(i)
        Pm1=zeros(ps,zs);
        Pw1=zeros(ps,zs);
        % kkk=kkk+1;
        clear Uv A Un a wi;
        Rxx=zeros(nel,nel);
        
        for j=1:nel
            for k=1:nel
                Rxx(j,k)=S(i,j,k);
            end
        end
        [Uv,A]=eig(Rxx);
        As=zeros(nel,nel);
        un=zeros(nel,nel);
        us=zeros(nel,nel);
       
        
        
        un(:,1:nel-M)=Uv(:,1:nel-M);
        
        us(:,1:M)=Uv(:,nel-M+1:nel);
        Un=un*un';
        sigma=0;
        for jj=1:nel-M
            sigma=sigma+1/(nel-M)*A(jj,jj);
        end
        
        vi=s(i);
        %wi=s1(i);
        wi=1;%ww(kkk);
        
        for p=1:ps
            
            for q=1:qs
                
                a=exp(-1i*2*pi*vi*tlib(:,p,q));
                Pm1(p,q)=(wi*(a'*a)/(a'*Un*a));
                Pw1(p,q)=(wi*(a'*Rxx*a)/(a'*a));
            end
        end
        %pause
        b=sort(reshape(Pm1,1,ps*qs));
        minp=mean(b(1:round(0.1*ps*qs)));
        maxp=max(max(Pm1));%old normalization
        Pm1=Pm1/max(max(Pm1));
        
        Pm=Pm+Pm1;
        Pw=Pw+Pw1;
        
    end
    
    
        