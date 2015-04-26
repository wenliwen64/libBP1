function matx=corrmat(x,sr,f)
[n m]=size(x);
matx=zeros(n,n);
NW=2;
  s=linspace(0,sr-sr/(m*sr),m*sr);
  fi=interp1(s,1:length(s),f);
for i=1:n
    for j=1:n
        [c ph ci ]=cmtm2(x(i,:),x(j,:),NW);
               
                 matx(i,j)=c(fi);
    end
end
                 