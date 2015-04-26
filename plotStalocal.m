function plotStalocal(DataFile)
load(DataFile);
figure(40);
close(40);
figure(40);
hold on;
r=ret.r;
hold on;
plot( r(:,1), r(:,2),'b.','MarkerSize',20);
xlabel('km');
ylabel('km');
axis equal
