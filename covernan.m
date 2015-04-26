function x=covernan(x)
% if there's Nan in a vector, replace the Nan by linear interpolation of
% the nearest points
a=[];
for i=1:length(x)
    if isnan(x(i))
        a=[a i];
    end
end
ind=1:length(x);
ind0=ind;
ind0(a)=[];
x0=x;
x0(a)=[];
x=interp1(ind0,x0,ind,'linear','extrap');