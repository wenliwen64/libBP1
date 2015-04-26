function ret=shiftAll(ret,shift)
  [n m]=size(ret.x);
 if length(shift)==1
     shift=ones(1,n)*shift;
 end
if isfield(ret,'x')
    for i=1:n
        ret.x(i,:)=specshift(ret.x(i,:),ret.sr*shift(i));
        
    end
end

if isfield(ret,'xori')
    for i=1:n
        ret.xori(i,:)=specshift(ret.xori(i,:),ret.sr*shift(i));
        
    end
end


if isfield(ret,'x1')
    for i=1:n
        ret.x1(i,:)=specshift(ret.x1(i,:),ret.sr*shift(i));
        
    end
end


% if isfield(ret,'recordtime')
% %     for i=1:n
%         ret.recordtime=ret.recordtime+shift/3600/24;
%         
% %     end
% end


% if isfield(ret,'t1')
% %     for i=1:n
%         ret.t1=ret.t1+shift;
%         
% %     end
% end

