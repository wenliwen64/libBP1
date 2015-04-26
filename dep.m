function  depth = dep(e,r,h)

if e==1
    depth=(r+1)/2*h(e);
else
    depth=sum(h(1:e-1))+(r+1)/2*h(e);
end
