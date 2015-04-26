function [xin yin]=peakfit2dlocal(Pm,minDist)
[mx my]=size(Pm);
[xin0 yin0]=localMaximum(Pm,minDist);
n=length(xin0);

for i=1:n
xin(i)=1;
yin(i)=1;
% find global maximum and extract 9-point neighbourship
% [v,p] = max(Z(:));
% [yp,xp]=ind2sub(sZ,p); 
yp=xin0(i)
xp=yin0(i)
if (yp==1)||(yp==mx)||(xp==1)||(xp==my)
    disp('Maximum position at matrix border. No subsample approximation possible.');
    P = [yp xp];
    return;
end
yp 
xp
K = Pm(yp-1:yp+1,xp-1:xp+1);

% approximate polynomial parameter
a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
%d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));

% (ys,xs) is subpixel shift of peak location relative to point (2,2)
ys = (6*b*c-8*a*f)/(16*e*a-9*b^2);
xs = (6*b*f-8*e*c)/(16*e*a-9*b^2);

P = [ys+yp xs+xp];
xin(i)=P(1);
yin(i)=P(2);
end