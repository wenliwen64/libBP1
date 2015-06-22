ep_lat = 28.1473; 
ep_lon = 84.7079;
bparea_span = [-0.5, 0.5; -0.5, 0.5];
bpgrid_x_vec = linspace(ep_lon + bparea_span(1, 1), ep_lon + bparea_span(1,2), 40);
bpgrid_y_vec = linspace(ep_lat + bparea_span(2, 1), ep_lat + bparea_span(2, 2), 40);
x_origin_10 = bpgrid_x_vec(10);
y_origin_10 = bpgrid_y_vec(10);
x_origin_30 = bpgrid_x_vec(30);
y_origin_30 = bpgrid_y_vec(30);
x_origin_20 = bpgrid_x_vec(20);
y_origin_20 = bpgrid_y_vec(20);
x_origin_ep = 0.5*(bpgrid_x_vec(20)+bpgrid_x_vec(21));
y_origin_ep = 0.5*(bpgrid_y_vec(20)+bpgrid_y_vec(21));
figure;
box on;
hold on;
plot(ep_lon, ep_lat, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'k');
plot(x_origin_10, y_origin_10, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
plot(x_origin_10, y_origin_30, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
plot(x_origin_30, y_origin_10, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
plot(x_origin_30, y_origin_30, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
plot(x_origin_20, y_origin_20, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
title('moho dip 15^o');
xlabel('lon(^o)');
ylabel('lat(^o)');

%bp_test_x = 84.72365; center
%bp_test_y = 28.1391;
%bp_test_x = 85.14576; %11 31

%bp_test_y = 27.58971;
% bp_test_x = 85.48077; +-0.5 resultc
% bp_test_y = 27.49618;

%bp_test_x = 85.22667;
%bp_test_y = 27.6239;
a_1010 = importdata('logfile_eutest_x10y10_dip15.dat');
[b, in] = max(a_1010(:, 4));
bp_test_y_1010 = a_1010(in, 2); % 28.60331;
bp_test_x_1010 = a_1010(in, 3); % 84.16469;
plot(bp_test_x_1010, bp_test_y_1010, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'm');

a_1030 = importdata('logfile_eutest_x10y30_dip15.dat');
[b, in] = max(a_1030(:, 4));
bp_test_y_1030 = a_1030(in, 2);
bp_test_x_1030 = a_1030(in, 3);
plot(bp_test_x_1030, bp_test_y_1030, 'd', 'MarkerSize', 20, 'MarkerFaceColor', 'm');

a_3010 = importdata('logfile_eutest_x30y10_dip15.dat');
[b, in] = max(a_3010(:, 4));
bp_test_y_3010 = a_3010(in, 2);
bp_test_x_3010 = a_3010(in, 3);
plot(bp_test_x_3010, bp_test_y_3010, 'd', 'MarkerSize', 20, 'MarkerFaceColor', 'm');

a_3030 = importdata('logfile_eutest_x30y30_dip15.dat');
[b, in] = max(a_3030(:, 4));
bp_test_y_3030 = a_3030(in, 2);
bp_test_x_3030 = a_3030(in, 3);
plot(bp_test_x_3030, bp_test_y_3030, 'd', 'MarkerSize', 20, 'MarkerFaceColor', 'm');

a_2020 = importdata('logfile_eutest_x20y20_dip15.dat');
[b, in] = max(a_2020(:, 4));
bp_test_y_2020 = a_2020(in, 2);
bp_test_x_2020 = a_2020(in, 3);
plot(bp_test_x_2020, bp_test_y_2020, 'd', 'MarkerSize', 20, 'MarkerFaceColor', 'm');

saveas(gcf, '/Users/lwen/Desktop/dip15_5locations', 'png');