close all;
N=256*8;
NW=30/2;
[E,V] = dpss(N,NW);
sr=40;
t=(1:N)/sr;
figure(20);
hold on;
for i=NW*2-2:NW*2
    plot(t,E(:,i),'color',rand(1,3));
end

figure(21);
 hold on;
for i=NW*2-2:NW*2
[Pxx,w] = pmtm(E(:,i),[]);
plot(w/(2*pi)*sr,Pxx,'color',rand(1,3));
 xlim([0 2])
end