% Ray tracing in a vertically heterogeneous medium with wave speed c(z)
% To shoot rays at given initial angles
% we solve the following ODE system
%   di/ds = dc/dz(z) * sin(i)
%   dx/ds = c(z)*sin(i)
%   dz/ds = c(z)*cos(i)
% where i = ray angle, counter-clockwise from Z down
%       s = travel-time
% The ODE is integrated by a midpoint method
%
% June 2008        J.P. Ampuero        ampuero@gps.caltech.edu

%I0S = [19.5:0.02:19.6];
I0S = [2:2:30];         % initial angles, in degrees
XSTA = 0; ZSTA = 0;        % initial point
%ZEND = 0;                % vertical position of final point
XEND = 5e3;                % horizontal position of final point
ZMAX = 20e3;                % stop ray tracing if the ray goes this deep
DS = 1e-2;                 % timestep

% example of velocity model c(z): an exponential distribution
CFUN = @(z) 2000 + (1-exp(-z/1000))*4000;
% its derivative with respect to z
CZFUN = @(z) exp(-z/1000)*4000/1000;

% for each initial angle, shoot a ray
for k=1:length(I0S),

 % initial conditions
  inc = I0S(k)*pi/180;
  x = XSTA;
  z = ZSTA;
  clear xray zray
  n=1;
  xray(n) =x;
  zray(n) =z;

 % Solve the ODE
 % Midpoint method to integrate Y'=F(Y) with step h:
 %  Y_pred = Y_n + h/2*F(Y_n)
 %  Y_(n+1) = Y_n + h*F(Y_pred)
 %
 % An alternative (also second order) is Heun's method:
 %  Y_pred = Y_n + h*F(Y_n)
 %  Y_(n+1) = Y_n + h*( F(Y_n) + F(Y_pred) )/2
 %
 % Or we could use ode45

  finish = 0;

  while ~finish,

   % Y_n
    i_old = inc;
    x_old = x;
    z_old = z;

   % F(Y_n)
    c = CFUN(z_old);
    cz = CZFUN(z_old);
    sini = sin(i_old);
    cosi = cos(i_old);
    is = cz*sini; 
    xs = c*sini; 
    zs = c*cosi; 

   % Y_pred = Y_n + h/2*F(Y_n)
    inc = i_old +0.5*DS*is;
    x = x_old +0.5*DS*xs;
    z = z_old +0.5*DS*zs;

   % F(Y_pred)
    c = CFUN(z);
    cz = CZFUN(z);
    sini = sin(inc);
    cosi = cos(inc);
    is = cz*sini; 
    xs = c*sini; 
    zs = c*cosi; 

   % Y_(n+1) = Y_n + h*F(Y_pred)
    inc = i_old +DS*is;
    x = x_old +DS*xs;
    z = z_old +DS*zs;

    n = n+1;
    xray(n)=x;
    zray(n)=z;

    %finish = z<ZEND | z>ZMAX;
    finish = x>XEND | z>ZMAX | z<ZSTA;
  end

  plot(xray/1e3,-zray/1e3,'-')
  hold all
end

hold off
axis equal
grid on
figure(5);
xxx=0:100:18000
plot(xxx,CFUN(xxx),'r',xxxx,ccc,'g');