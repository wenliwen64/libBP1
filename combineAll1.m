function ret=orderAll(ret1,ret2)


if isfield(ret1,'x')&&isfield(ret2,'x')
ret.x=[ ret1.x; ret2.x];
end
if isfield(ret1,'xori')&&isfield(ret2,'xori')
ret.xori=[ret1.xori ;ret2.xori];
end

if isfield(ret1,'r')&&isfield(ret2,'r')
ret.r=[ret1.r ;ret2.r];
end
if isfield(ret1,'nm')&&isfield(ret2,'nm')
ret.nm=[ret1.nm ;ret2.nm];
end

if isfield(ret1,'lat')&&isfield(ret2,'lat')
ret.lat=[ret1.lat ret2.lat];
end
if isfield(ret1,'lon')&&isfield(ret2,'lon')
ret.lon=[ret1.lon ret2.lon];
end
if isfield(ret1,'rdis')&&isfield(ret2,'rdis')
ret.rdis=[ret1.rdis ret2.rdis];
end
if isfield(ret1,'az')&&isfield(ret2,'az')
ret.az=[ret1.az ret2.az];
end
if isfield(ret1,'time')&&isfield(ret2,'time')
ret.time=[ret1.time ret2.time];
end
if isfield(ret1,'t1')&&isfield(ret2,'t1')
ret.t1=[ret1.t1 ret2.t1];
end
if isfield(ret1,'x1')&&isfield(ret2,'x1')
    ret.x1=[ret.x1; ret.x2];
end
if isfield(ret1,'timeshift')&&isfield(ret2,'timeshift')
    ret.timeshift=[ret1.timeshift ret2.timeshift];
end


    
