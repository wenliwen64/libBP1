figure(1);
%M=1;
nel=12;
kv=@(f0,ux,uy) 2*pi*f0*[ux uy (1/2^2-ux^2-uy^2)^0.5];% wave vector
v=@(r,k) exp(-1i*(r*k'));%steering vector
r=zeros(nel,3);%coordinates of the array element
r(:,1)=[-371 -298 -192 -23 0 23 177 213 260 98 -72 -131 ]/1000;
r(:,2)=[-302 -208 -284 -13 0 -9 98 235 411 213 328 393 ]/1000;
r(:,3)=[-25 -24 -16  -4 0 2 18 12 3 -1 -16 -4 ]/1000;

ps=100;
qs=100;
Xm=zeros(ps,qs);
Ym=zeros(ps,qs);
Pm=zeros(ps,qs);
ux=linspace(-1,1,ps);
uy=linspace(-1,1,qs);

for q=1:qs
    Xm(:,q)=ux;
end
for p=1:ps
    Ym(p,:)=uy';
end

   for p=1:ps
        for q=1:qs
            
            a=v(r,kv(10,ux(p),uy(q)));
            Pm(p,q)=(sum(a))^2;
        end
   end
    
   
   contourf(Xm,Ym,Pm,15);