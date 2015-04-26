figure(1);

red=[ 1 0 0];
green=[0 1 0];
blue=[0 0 1];

r=interp1([1 10 20],[1 0 0],1:20)
g=interp1([1 10 20],[0 1 0],1:20)
b=interp1([1 10 20],[0 0 1],1:20)

for i=1:20
    hold on;
    plot(i,i,'color',[r(i) g(i) b(i)],);
end