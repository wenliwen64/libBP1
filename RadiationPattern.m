% RADIATIONPATTERN radiation pattern of a double couple point source
%
% SYNTAX        U = RadiationPattern(str,dip,rak, AIN, AZM)
%
% INPUTS        str        strike angle
%                dip        dip angle
%                rak        rake angle
%                AIN        take-off angle
%                AZM        azimuth angle
%
% NOTE                All angles are in degrees. 
%                Definitions follow usual seismological conventions.
%                The inputs can be vectors of same length N.
%
% OUTPUTS        U        radiation pattern, size = [N,3]
%                        U(:,1) = P
%                        U(:,2) = SV 
%                        U(:,3) = SH
%
function U = RadiationPattern(A,B,C,D,E,R)

rad = pi/180;
A = A(:)*rad;        % strike (phi_f)         
B = B(:)*rad;        % dip (delta)
C = C(:)*rad;         % rake (lambda) 
D = D(:)*rad;        % AIN
E = E(:)*rad;        % AZM

U = zeros(max([length(A),length(B),length(C),length(D),length(E)]),3);

% P
U(:,1) = cos(C).*sin(B).*sin(D).^2.*sin(2*(E-A)) ...
    - cos(C).*cos(B).*sin(2*D).*cos(E-A) ...
    + sin(C).*sin(2*B).*(cos(D).^2-sin(D).^2.*sin(E-A).^2) ...
    + sin(C).*cos(2*B).*sin(2*D).*sin(E-A);

% SV
U(:,2) = sin(C).*cos(2*B).*cos(2*D).*sin(E-A) ...
    - cos(C).*cos(B).*cos(2*D).*cos(E-A) ...
    + .5*cos(C).*sin(B).*sin(2*D).*sin(2*(E-A)) ...
    - .5*sin(C).*sin(2*B).*sin(2*D).*(1 + sin(E-A).^2);

% SH
U(:,3) = cos(C).*cos(B).*cos(D).*sin(E-A) ...
    + cos(C).*sin(B).*sin(D).*cos(2*(E-A)) ...
    + sin(C).*cos(2*B).*cos(D).*cos(E-A) ...
    - .5*sin(C).*sin(2*B).*sin(D).*sin(2*(E-A));

if exist('R','var')
  R = R(:);
  U = U./[R R R];
end