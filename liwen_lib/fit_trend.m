close all;
trend_y = y(3501:4000);
trend_func = @(x,xdata)x(1)*exp(-xdata*x(2))+x(3);
x0 = [4.8, 0.1, 3.5];
t = 1:500;
[x, resnorm, ~, exitflag, output]=lsqcurvefit(trend_func, x0,t,trend_y)
figure;
hold on;
plot(t, trend_y);
plot(t, trend_func(x, t));

figure;
hold on;
plot(1:1000, y(3501:4500));
plot(1:1000, trend_func(x, 1:1000));


trend_y2 = y(3501:4500);
% for i = 1:numel(trend_y2)
    trend_y2 = trend_y2 - trend_func(x, 1:1000); 
% end
figure;
plot(trend_y2);