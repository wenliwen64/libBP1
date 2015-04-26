%
% plot_basic.m
% Yoshihiro Kaneko, 16-Oct-2007
% 
% This script is used for matlab tutorial
%

clear all       % clear all variables in memory
close all       % close all open figures

% generate 50 linearly-spaced points, ranging from -2 to 2
x = linspace(-2,2,50);    

% generate "data points" using sine function
data = sin(2*pi*x);

% plot continuous curve
figure;
subplot(3,2,1);            % make 3 by 2 sub-panels
plot(x,data)               
xlabel('x');
ylabel('data');
title('sine function');

orient tall,wysiwyg;      % make the size of the figure on the screen to equal
                           % the size of the figure that would be printed 

%%% same example with additional plotting commands %%%

subplot(3,2,2);            
plot(x,data,'k.-')                        
xlabel('x','fontsize',12);
ylabel('data','fontsize',12);
set(gca,'xtick',[-1 -0.3 0 0.8 1.2]); % control tick spacing for x-axis
set(gca,'ytick',[-1:0.25:1]);         % control tick spacing for y-axis
axis([-1.5 1.5 -1 1]);                % control axis scaling and appearance


%%% another example: log scale is used for the Y-axis %%%

sigma_h = exp(exp(x));

% plot discrete data points
subplot(3,2,3);
semilogy(x,sigma_h,'m.')
xlabel('x','fontname','Times','fontsize',16);
ylabel('\sigma_h','fontname','Times','fontsize',16);
title('e^{e^{x}}');
set(gca,'ydir','reverse');



%%% another example: Log-log scale plot %%%

% generate 20 logrithmically-spaced points, ranging from 10^-10 to 10^10
x = logspace(-10,10,20);
func1 = 4*x;
func2 = 900*x;

% plot discrete data points
subplot(3,2,4);
loglog(x,func1,'r-.',x,func2,'b--')
xlabel('x');
ylabel('y','fontname','Times','fontsize',16);
legend('{\it f}_1','{\it f}_2','location','southeast');


%%% another example: 3-D colored surface plot  %%%

clear all;
subplot(3,2,[5 6]);
load topo.txt    % load topographys data called "topo.txt"
d=topo;

surfc(topo);
xlabel('West-East Distance (m)');
ylabel('North-South Distance (m)');
zlabel('Elevation (m)');
title('Topography');
colorbar;


