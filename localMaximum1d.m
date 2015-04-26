function [x kk] = localMaximum1d(Pm)
m=length(Pm);
x=0;

k=1;
for i=2:m-1
  
        if Pm(i)>max([Pm(i-1) Pm(i+1)  ])
            x(k)=i;
           
            k=k+1;
        end
 
end
kk=k-1;