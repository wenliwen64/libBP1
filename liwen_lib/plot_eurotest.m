ep_lat = 28.1473; 
ep_lon = 84.7079;
bparea_span = [-0.5, 0.5; -0.5, 0.5];
bpgrid_x_vec = linspace(ep_lon + bparea_span(1, 1), ep_lon + bparea_span(1,2), 40);
bpgrid_y_vec = linspace(ep_lat + bparea_span(2, 1), ep_lat + bparea_span(2, 2), 40);
x_origin = bpgrid_x_vec(15);
y_origin = bpgrid_y_vec(25);
figure;
hold on;
plot(ret.lon0, ret.lat0, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'k');
plot(x_origin, y_origin, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
%bp_test_x = 84.72365; center
%bp_test_y = 28.1391;
%bp_test_x = 85.14576; %11 31

%bp_test_y = 27.58971;
% bp_test_x = 85.48077; +-0.5 resultc
% bp_test_y = 27.49618;

%bp_test_x = 85.22667;
%bp_test_y = 27.6239;
a = importdata('logfile_eutest_x15y25_dip20.dat');
[b, in] = max(a(:, 4));
bp_test_y = a(in, 2) % 28.60331;
bp_test_x = a(in, 3) % 84.16469;
plot(bp_test_x, bp_test_y, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'm');