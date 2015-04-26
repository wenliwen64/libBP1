% LOVE_ANALYTIC computes the exact solution for Homework #1 (Love modes)
% 
% SYNTAX        [u,k] = love_analytic(z,m)
%
% INPUTS        z        depth values (positive, can be a vector or a matrix)
%                m        mode number, 0=fundamental or 1=first overtone
%
% OUTPUTS        u         displacement amplitude at requested depths
%                        normalized such that u(z=0)=1
%                k        horizontal wavenumber
% 
function [u,k] = love_analytic(z,m)

% model parameters:
RHO = 2000;
H1 = 400;
H2 = 3600;
C1 = 1000;
C2 = 2000;
FREQ = 2.5;

w = 2*pi*FREQ;
MU1 = RHO*C1^2;
MU2 = RHO*C2^2;

% Find one of the two largest eigenvalues (k0 or k1)
% The dispersion relation is
%   tan(k1*H1)/tanh(g1*H2) = mu2*g2/(mu1*k1)
% where
%   k1 = sqrt(w^2/C1^2-k^2);
%   g2 = sqrt(k^2-w^2/C2^2);
%
% Note: H2 is so large that tanh(g1*H2) is very close to 1 and
% the problem almost reduces to a half-space with a shallow layer 
% (eq. 7.6 of Aki & Richards, 2002)
% The range to search for k for each mode
% is based on Figure 7.2 of Aki and Richards (2002)
%
if m<0 | m>1, error('Only m=0 or 1 are allowed'); end
krange = w/C1*sqrt(1-([m+1/2 m]*pi/w*C1/H1).^2) ;
[k,fval,exitflag] = fzero(@fun,krange);

% evaluate the slip profile
k1 = sqrt(w^2/C1^2-k^2);
g2 = sqrt(k^2-w^2/C2^2);
u = (z<H1) .* cos(k1*z)  ...
  + (z>=H1) .* sinh(g2*(H1+H2-z))* cos(k1*H1)/sinh(g2*H2);

%-------
function f=fun(kx)
k1 = sqrt(w^2/C1^2-kx^2);
g2 = sqrt(kx^2-w^2/C2^2);
f = tan(k1*H1)*MU1*k1 -tanh(g2*H2)*MU2*g2;
end

end
