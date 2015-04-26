function [Y,I] = csort(X,method,dim,mode);
% CSORT   Sort complex 2D-arrays.
%     For vectors, SORT(X) sorts the elements of X in ascending order.
%     For matrices, SORT(X) sorts each column of X in ascending order.
%     The standard method is a lexical sort, first sort by real part,
%     then by imaginary part.
% 
%     Y = CSORT(X,METHOD,DIM,MODE)
%     has three optional parameters.  
%     METHOD selects the method of the sort
%        'lexi' lexical sort, sort first by REAL(X), then by IMAG(X)
%        'angle' sort by phase angle, ANGLE(X)
%        'abs' sort by magnitude, ABS(X)
%        'real' sort by real part, REAL(X)
%        'imag' sort by imaginary part, imag(X)
%        'absangle' sort first by ABS(X), then ANGLE(X) (Matlab standard)
%     DIM selects a dimension along which to sort.
%     MODE selects the direction of the sort
%        'ascend' results in ascending order
%        'descend' results in descending order
%     The MODE option is not valid for lexical sorting.
%     The result is in Y which has the same shape and type as X
% 
%     [Y,I] = CSORT(X,METHOD,DIM,MODE) also returns an index matrix I.
%     If X is a vector, then Y = X(I).  
%     If X is an m-by-n matrix and DIM=1, then
%         for j = 1:n, Y(:,j) = X(I(:,j),j); end



%  Author: Peter Bodin <pbodin@kth.se>
%  2006-01-08 
error(nargchk(1,4,nargin,'struct'));
if nargin < 4 | isempty(mode), mode = 'ascend'; end
if nargin < 3 | isempty(dim), dim = 1; end
if nargin < 2 | isempty(method), method = 'lexi'; end

switch method
    case 'lexi'             % lexical sort, first by real part then imaginary part
        if dim ~= 1
           [Y,I] = csort(X',method,1);
           Y = Y';
           return
        end
        [I,I] = sortrows([real(X) imag(X)],[1 2]);
    case 'angle'            % sort by phase angle
        [I,I] = sort(imag(log(X)),dim,mode);
    case 'real'             % sort by real part only
        [I,I] = sort(real(X),dim,mode);
    case 'imag'             % sort by imaginary part only
        [I,I] = sort(imag(X),dim,mode);
    case 'abs'             % sort by magnitude
        [I,I] = sort(abs(X),dim,mode);
    case 'absangle'        % sort first by magnitude then angle angle (Matlab standard)
        [I,I] = sort(X,dim,mode);
    case 'otherwise'
        error('CSORT: Allowed methods are: lexi, angle, real, imag, abs or magnphase')
end

Y = X(I,:);