%
% function x = randomvec(xmin, xmax, n)
% Carl Tape, 13-June-2005
% printed xxx
%
% This generates a vector of n random numbers in the range [xmin,xmax].
%

function x = randomvec(xmin0, xmax0, n)

xmin = min([xmin0 xmax0]);
xmax = max([xmin0 xmax0]);

x = (xmax - xmin)*rand(n,1) + xmin;

%================================================
