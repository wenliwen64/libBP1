%used in plotEEWweb
function [x z]=quadratic3d(y,x0,x1,y1,z0,z01)
L2=(x1-x0)^2+y1^2;
x=x0+y/y1*(x1-x0);
z=(z01-z0)*((x-x0).^2+y.^2)/L2+z0;
end