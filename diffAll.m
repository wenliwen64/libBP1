function ret=diffAll(ret)
[n m]=size(ret.x);
xi=zeros(n,m);
for i=1:n
    xi(i,1:end-1)=diff(ret.x(i,:));
    xi(i,end)=xi(i,end-1);
end
ret.x=xi;