function ret=densityWeight(ret,radius)
n=ret.n;
weight=zeros(1,n);
for i=1:n
    count=0;
    for j=1:n
        if distance(ret.lat(i),ret.lon(i),ret.lat(j),ret.lon(j))<radius
            count=count+1;
        end
        
    end
    weight(i)=1/count;
    ret.x(i,:)=weight(i)*ret.x(i,:);
end
            ret.weight=weight;