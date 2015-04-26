clear all;
close all;

kv=@(f0,ux,uy) 2*pi*f0*[ux uy];% wave vector
v=@(r,k) exp(1i*(r*k'));%steering vector

f=10;
nel=400;
r=zeros(nel,2);
r(:,1)=linspace(-3,3,nel);
sm=0.2;




    ps=100;
    qs=100;
    Xm=zeros(ps,qs);
    Ym=zeros(ps,qs);
    Pm=zeros(ps,qs);
    ux=linspace(-sm,sm,ps);
    uy=linspace(-sm,sm,qs);
    
    for q=1:qs
        Xm(:,q)=ux;
    end
    for p=1:ps
        Ym(p,:)=uy';
    end
    
    for p=1:ps
        for q=1:qs
            
            a=v(r,kv(f,ux(p),uy(q)));
            Pm(p,q)=((abs(sum(a)))^2);
        end
    end
 figure(1);
 subplot(121);
 pcolor(ux,uy,Pm);
 shading flat;
 colorbar;
subplot(122);
plot(ux,Pm(:,1));







