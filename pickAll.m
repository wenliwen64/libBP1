function [ret ]=pickAll(ret,tb,te)
[n m]=size(ret.x);
t3=(0:length(ret.x)-1)/ret.sr;
h=figure(95);
ret.tp=zeros(1,n);
i=1;
while(i<=n)
    display(num2str(i))
plot(t3,ret.x(i,:),[ret.tp(i) ret.tp(i)],[-max(ret.x(i,tb*ret.sr:te*ret.sr)) max(ret.x(i,tb*ret.sr:te*ret.sr))],'r');
title(['N' num2str(i) ' ' ret.nm(i,:) ]);
xlim([tb te]);
[px, py, button]=ginput(1);
    display(num2str(button));

if button==28&&i>=2
      display(num2str(button));
    i=i-1;
    continue;
elseif button==29&&i<=n
      display(num2str(button));
    i=i+1;
    continue;
elseif button==1
     ret.tp(i)=px;
           display(['i3 ' num2str(i)]);

     i=i+1;
           display(['i4 ' num2str(i)]);

else
title(['wrong input']);    
end

% ret.tp(i)=px;
end
for i=1:n
    ret.x(i,:)=specshift(ret.x(i,:),(ret.tp(i)-10)*ret.sr)';
    ret.xori(i,:)=specshift(ret.xori(i,:),(ret.tp(i)-10)*ret.sr)';

end

