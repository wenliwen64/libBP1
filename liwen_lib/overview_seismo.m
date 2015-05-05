%load mainshock.v2.mat
%close all;
%ret = new_str;
nsta = ret.n;
figure(2);
hold all;
for i = 1:400
    %plot(ret.xori(i, 20000:35000)/std(ret.xori(i, 1000:4000)) + i*100);
        plot(ret.xori(i, 1:20000)/std(ret.xori(i, 1000:4000)) + i*100);

end