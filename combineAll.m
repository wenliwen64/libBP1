function ret=orderAll(ret,ret1)
  

if isfield(ret,'x')
N=min([length(ret.x) length(ret1.x)]);% using the shorter length of x
ret.x=[ret.x(:,1:N); ret1.x(:,1:N)];
end

if isfield(ret,'xori')
N=min([length(ret.x) length(ret1.x)]);
ret.xori=[ret.xori(:,1:N); ret1.xori(:,1:N)];
end

if isfield(ret,'r')
ret.r=[ret.r; ret1.r];
end
if isfield(ret,'nm')
ret.nm=[ret.nm ;ret1.nm];
end

if isfield(ret,'lat')
ret.lat=[ret.lat ret1.lat];
end
if isfield(ret,'lon')
ret.lon=[ret.lon ret1.lon];
end
if isfield(ret,'rdis')
ret.rdis=[ret.rdis ret1.rdis];
end
if isfield(ret,'az')
ret.az=[ret.az ret1.az];
end
% if isfield(ret,'time')
% ret.time=[ret.time ;ret1.time];
% end
if isfield(ret,'t1')
ret.t1=[ret.t1 ret1.t1];
end

if isfield(ret,'recordtime')
ret.recordtime=[ret.recordtime ret1.recordtime];
end
if isfield(ret,'x1')
N=min([length(ret.x1) length(ret1.x1)]);
ret.x1=[ret.x1(:,1:N); ret1.x1(:,1:N)];
end
if isfield(ret,'timeshift1')
    ret.timeshift1=[ret.timeshift1; ret1.timeshift1];
end

[ret.n ret.m]=size(ret.x);
    
