function ret=intAll(ret)
[n m]=size(ret.x);
xi=zeros(n,m);
for i=1:n
    xi(i,:)=cumsum(ret.x(i,:));
end
ret.x=xi;