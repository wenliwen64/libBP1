clear;
load wave;
lat=wave(:,1);
lon=wave(:,2);
time=wave(:,3);
anomaly=wave(:,4);
residue=wave(:,5);
totalfield=wave(:,7);
figure;
count=817;
for i=1:1:count
    if residue(i)>9000
    anomaly(i)=0;
    residue(i)=0;
    end
end

lati=linspace(-5,0.0002*817-5,817);

plot(lati,totalfield,'r');
