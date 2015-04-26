function [retc retc1]=matchAll(ret,ret1)
    [lon B B1]=intersect(ret.lon,ret1.lon,'sorted');
    retc=orderAll(ret,0,B);
    retc1=orderAll(ret1,0,B1);
end
    

    
