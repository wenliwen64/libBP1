function ret=decimateAll(ret,sr)
[n m]=size(ret.x);
N=floor(ret.sr/sr);
if N<2
    return
end
m1=length(decimate(ret.x(1,:),N));
x1=zeros(n,m1);
x2=zeros(n,m1);
for i=1:n
    x1(i,:)=decimate(ret.x(i,:),N);
    x2(i,:)=decimate(ret.xori(i,:),N);
end
ret.x=x1;
ret.xori=x2;
ret.sr=ret.sr/N;