%
% test_getsubset.m
% Carl Tape, 10-Oct-2006
% modified by Yoshihiro Kaneko, 17-Oct-2007
%
% This program shows how to extract a sub-box of datapoint indices
% from a set of (irregular) xy data.
%
% Note that for plotting ANY function, Matlab requires regular
% spacing between gridpoints. Thus, here we use the griddata.m
% function for plotting purposes.
%
% calls getsubset.m, gridvec.m, randomvec.m, griddataXB.m
% called by xxx
%

clear           % clear all variables in memory
close all       % close all open figures
format short    % limit command-window output to single precision
format compact  % reduce spacing in command-window output
%clc            % erase the command window

%------------------------------
% generate the x-y test points (this can be done in several ways)

% designate the bounds of region
ax0 = [10 100 20 200];
xmin = ax0(1); xmax = ax0(2);
ymin = ax0(3); ymax = ax0(4);

numx = 50;      % determines spacing in the horizontal dimension

% random gridpoints
xvec = randomvec(xmin,xmax,
    numx*numx);
yvec = randomvec(ymin,ymax,numx*numx);

num = length(xvec);

% assign z = f(x,y) values, where f is a test function
w = 2*pi/40;
zvec = sin(xvec*w) .* sin(yvec*w);

%------------------------------

% compute a regular grid for plotting
npts = 100;
[X,Y,Z] = griddataXB(xvec,yvec,zvec,npts,'cubic');

figure; nr=2; nc=1;
clim = [-1 1];

subplot(nr,nc,1); hold on;
pcolor(X,Y,Z); shading interp;
plot(xvec,yvec,'k+','markersize',2);
xlabel(' x'); ylabel(' y');

fontsize(12), orient landscape, wysiwyg;    % set fontsize
                                            % set paper orientation to be "landscape"
                                                
%------------------------------
% select subset 'window'

xpmin = 14.4; xpmax = 48.0;
ypmin = 42.0; ypmax = 83.2;

% bounding box for subset domain
xbox = [xpmin xpmax xpmax xpmin xpmin];
ybox = [ypmin ypmin ypmax ypmax ypmin];
plot(xbox, ybox, 'r', 'linewidth', 3);

% KEY COMMAND: THIS GETS THE DESIRED INDICES
[ix,x2,y2] = getsubset(xvec,yvec,[xpmin xpmax ypmin ypmax]);
nump = length(ix);

plot(x2,y2,'k.'); caxis(clim);
xlabel(' x'); ylabel(' y');
title(['test function plotted on irregular datapoints : ' num2str(num) ' total points']);

%------------------------------
% plot the function sampled by the subset of datapoints

npts = 40;
[Xr,Yr,Zr] = griddataXB(x2,y2,zvec(ix),npts,'cubic');

subplot(nr,nc,2); hold on;
pcolor(Xr,Yr,Zr); shading interp; caxis(clim);
plot(x2,y2,'k.');
plot(xbox, ybox, 'r', 'linewidth', 3);
xlabel(' x'); ylabel(' y');
title(['subset of points extracted from above (' num2str(nump) ' / ' num2str(num) ')']);

fontsize(12), orient landscape, wysiwyg;    % set fontsize
                                            % set paper orientation to be "landscape"

%=====================================================
