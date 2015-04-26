function ret=removeAll(ret,re_i)
 a=[];
[n m]=size(ret.x);
 for i=1:n
     if ~ismember(i,re_i)
         a=[a i];
     end
 end
B=a;
ret.n=length(B);
if isfield(ret,'x')
ret.x=ret.x(B,:);
end

if isfield(ret,'xori')
ret.xori=ret.xori(B,:);
end

% if isfield(ret,'x1')
%     B
%     size(ret.x1)
% ret.x1=ret.x1(B,:);
% end
if isfield(ret,'r')
ret.r=ret.r(B,:);
end
if isfield(ret,'nm')
ret.nm=ret.nm(B,:);
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
if isfield(ret,'time')
ret.time=ret.time(B,:);
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
