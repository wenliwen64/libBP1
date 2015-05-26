% this script plot test points for Australia BP test

% Original data input:
close all;
ep_lat = 28.1473; 
ep_lon = 84.7079;
bparea_span = [-1, 1; -1, 1];
bpgrid_x_vec = linspace(ep_lon + bparea_span(1, 1), ep_lon + bparea_span(1,2), 40);
bpgrid_y_vec = linspace(ep_lat + bparea_span(2, 1), ep_lat + bparea_span(2, 2), 40);
test_point_ind_x = [10, 10, 15, 25, 30, 30];
test_point_ind_y = [10, 30, 25, 15, 10, 30];

test_point_x = bpgrid_x_vec(test_point_ind_x);
test_point_y = bpgrid_y_vec(test_point_ind_y);

figure;
plot(test_point_x, test_point_y, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');

% Bp data point
for i = 1:numel(test_point_x)
    filename = ['logfile_eutest_' num2str(test_point_ind_x(i)) '_' num2str(test_point_ind_y(i)) '.dat'];
    filename
    a = importdata(filename);
    [b, ind] = max(a(:, 4));
    bp_point_x(i) = a(ind, 3);
    bp_point_y(i) = a(ind, 2);
end
hold on;
plot(bp_point_x, bp_point_y, 'd', 'MarkerSize', 20,  'MarkerFaceColor', 'm');   

plot(ep_lon, ep_lat, 'd', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% a = importdata('logfile_autest.dat');
% [b, ind] = max(a(:, 4));
% plot(a(ind, 3), a(ind, 2), 'd', 'MarkerSize', 10 ,'MarkerFaceColor', 'm');
axis([ep_lon-1 ep_lon+1 ep_lat-1 ep_lat+1]);
xlabel('lon');
ylabel('lat');
title('Europe Test');
tand(155);
x_moho_line = ep_lon-1:0.02:ep_lon+1;
y_moho_line = tand(155)*(x_moho_line - ep_lon) + ep_lat;
plot(x_moho_line, y_moho_line, '*', 'MarkerFaceColor', 'green');