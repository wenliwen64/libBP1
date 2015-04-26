function [x y] = localMaximum1(Pm)
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