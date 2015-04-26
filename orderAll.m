function ret=orderAll(ret,dist,ind)
if dist==0
    B=ind;
    
else
    if nargin<3% no opr speciified, use default
        ind=1:length(dist);
    end
    
    tmp=sortrows([ind' dist(ind)'],2);
    
    B=tmp(:,1);
end
ret.n=length(B);
size(B)
B(end)
if isfield(ret,'x')
ret.x=ret.x(B,:);
end
if isfield(ret,'xori')
ret.xori=ret.xori(B,:);
end

if isfield(ret,'r')
ret.r=ret.r(B,:);
end
if isfield(ret,'nm')
ret.nm=ret.nm(B,:);
end
if isfield(ret,'mf')
ret.mf=ret.mf(:,B);
end
if isfield(ret,'wnorm')
ret.wnorm=ret.wnorm(B);
end
if isfield(ret,'lat')
ret.lat=ret.lat(B);
end
if isfield(ret,'lon')
ret.lon=ret.lon(B);
end
if isfield(ret,'rdis')
ret.rdis=ret.rdis(B);
end
if isfield(ret,'az')
ret.az=ret.az(B);
end
% if isfield(ret,'time')
% ret.time=ret.time(B,:);
% end
if isfield(ret,'recordtime')
ret.recordtime=ret.recordtime(B);
end
if isfield(ret,'t1')
ret.t1=ret.t1(B);
end
if isfield(ret,'x1')
    ret.x1=ret.x1(B,:);
end
if isfield(ret,'timeshift')
    ret.timeshift=ret.timeshift(B);
end
if isfield(ret,'timeshiftall')
    ret.timeshiftall=ret.timeshiftall(B);
end
if isfield(ret,'nf')
    ret.nf=ret.nf(B);
end

    
