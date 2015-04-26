function [x y z] = localMaximumfit(Pm)
[m n]=size(Pm);
Pm1=zeros(m+2,n+2);
Pm1(2:m+1,2:n+1)=Pm;
x=0;
y=0;
k=1;
for i=2:m+1
    for j=2:n+1
         if Pm1(i,j)>max([Pm1(i-1,j-1) Pm1(i,j-1) Pm1(i-1,j) Pm1(i+1,j)  Pm1(i,j+1) Pm1(i+1,j+1)])
            x(k)=i-1;
            y(k)=j-1;
            k=k+1;
        end
    end
end

x=x';
y=y';

nn=length(x);
for i=1:nn
    yp=x(i)+1
    xp=y(i)+1
K = Pm1(yp-1:yp+1,xp-1:xp+1);

% approximate polynomial parameter
a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
%d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));

% (ys,xs) is subpixel shift of peak location relative to point (2,2)
ys = (6*b*c-8*a*f)/(16*e*a-9*b^2);
xs = (6*b*f-8*e*c)/(16*e*a-9*b^2);

P = [ys+yp xs+xp];
x(i)=real(P(1)-1);
y(i)=real(P(2)-1);
z(i)=real(interp2(1:m,1:n,Pm,x(i),y(i)));
end

tmp=sortrows([x y z'],-3);
x=tmp(:,1);
y=tmp(:,2);
z=tmp(:,3);