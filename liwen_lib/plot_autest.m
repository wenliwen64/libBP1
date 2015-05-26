%x_origin = bpgrid_x_vec(25); % ret.lon0-1+2/40*20;
%y_origin = bpgrid_y_vec(15); % ret.lat0-1+2/40*20;
figure;
hold on;
plot(ret.lon0, ret.lat0, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'k');
%plot(x_origin, y_origin, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
%bp_test_x = 84.72365; center
%bp_test_y = 28.1391;
%bp_test_x = 85.14576; %11 31

%bp_test_y = 27.58971;
% bp_test_x = 85.48077; +-0.5 resultc
% bp_test_y = 27.49618;

%bp_test_x = 85.22667;
%bp_test_y = 27.6239;
a = importdata('logfile_autest.dat');
[b, in] = max(a(:, 4));
bp_test_y = a(in, 2) % 28.60331;
bp_test_x = a(in, 3) % 84.16469;
plot(bp_test_x, bp_test_y, 'd', 'MarkerSize', 20, 'MarkerFaceColor', 'm');