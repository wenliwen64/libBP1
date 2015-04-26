eqinfo1 = importdata('../eqinfo.dat');
eqinfo2 = importdata('./eqinfo_2014.dat');
begint = datenum([2014 4 1 23 0 0]);
endt = datenum([2014 4 2 23 0 0]);
figure(1);
hold on;
I = find(eqinfo1(:, 7) >= begint & eqinfo1(:,7) < endt);
II = find(eqinfo2(:, 7) >= begint & eqinfo2(:,7) < endt);

plot(eqinfo1(I,9), eqinfo1(I,8), 'b.', 'Markersize', 20);
plot(eqinfo2(II,9), eqinfo2(II,8), 'r*', 'Markersize', 20);