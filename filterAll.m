%function ret=filterAll(ret,fl,fh)
function ret=filterAll(ret,fl,fh)
ret.x=ret.xori;
[n m]=size(ret.x);
[BB,AA]=butter(4,[fl fh]/(ret.sr/2));
for i=1:n
    
    ret.x(i,:)=filter(BB,AA,ret.x(i,:));
end