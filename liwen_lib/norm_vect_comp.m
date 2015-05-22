% this script used to find the normal vector of a cern plane from strike
% and dip
strike_angle = 295;
dip_angle = 11;
x0 = 1*cosd(strike_angle);
y0 = 1*sind(strike_angle);
z0 = 0;

x1 = 1*cosd(dip_angle)*cosd(strike_angle+90);
y1 = 1*cosd(dip_angle)*sind(strike_angle+90);
z1 = -1*sind(dip_angle);

a = [x0, y0, z0];
b = [x1, y1, z1];
normal_vec = cross(a,b);