x_origin = ret.lon0+0.5;
y_origin = ret.lat0-0.5;
figure;
hold on;
plot(ret.lon0, ret.lat0, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'k');
plot(x_origin, y_origin, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'b');
%bp_test_x = 84.72365; center
%bp_test_y = 28.1391;
%bp_test_x = 84.1948; %11 31

%bp_test_y = 28.6525;
% bp_test_x = 85.48077; +-0.5 result
% bp_test_y = 27.49618;

bp_test_x = 85.22667;
bp_test_y = 27.6239;
plot(bp_test_x, bp_test_y, 'd', 'MarkerSize', 20,'MarkerFaceColor', 'm');