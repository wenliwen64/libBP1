function [ret ret1]=matchshiftAll(ret,ret1)
    [lon,I,J]=unique(ret.lon);
    ret=orderAll(ret,0,I);
    [lon,I1,J]=unique(ret1.lon);
    ret1=orderAll(ret1,0,I1);
    
    
    [lon B B1]=intersect(ret.lon,ret1.lon);
    
    ret=orderAll(ret,0,B);
    ret1=orderAll(ret1,0,B1);
    nel=length(ret.lon);
    if isfield(ret,'x')
        for i=1:nel
            ret.x(i,:)=specshift(ret.x(i,:),(ret1.timeshiftall(i)-ret.timeshiftall(i))*ret.sr)/ret1.nf(i);
        end
    end
    if isfield(ret,'xori')
        for i=1:nel
            ret.xori(i,:)=specshift(ret.xori(i,:),(ret1.timeshiftall(i)-ret.timeshiftall(i))*ret.sr)/ret1.nf(i);
        end
    end
    ret.timeshiftall=ret1.timeshiftall;
    ret=orderAll(ret,ret.rdis,1:nel);
    ret1=orderAll(ret1,ret1.rdis,1:nel);

end
    

    
