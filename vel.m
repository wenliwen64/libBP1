function  velc = vel(e,r,h)

if e==1
    z=(r+1)/2*h(e);
else
    z=sum(h(1:e-1))+(r+1)/2*h(e);
end
if z<400
    velc=1000;
else
    velc=2000;
end

    
