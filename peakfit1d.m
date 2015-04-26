function [val ind] = peakfit1d(x) 
%find the quandratic fitting peak of the vector x
[val ind]=max(x);
if ind==1||ind==length(x)
    return
end
[d_ind,val,a] = qint(x(ind-1),x(ind),x(ind+1));
ind=ind+d_ind;