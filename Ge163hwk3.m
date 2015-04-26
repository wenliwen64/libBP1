pe=5.5153e3
pw=1000
pc=2700
pm=3300
hl=1000
t=100e3
Re=6370e3
dN=(1:99)*0;
for l=2:100
    dN(l-1)=-3/pe*(1/(2*l+1)*hl*(pc-pw)*(1-(1-t/Re)^(l+2)))
end
plot(dN);
title('dN_l vs l ( t = 100 km )');
xlabel('degree l');
ylabel('geoid anomaly dN_l')
N=sum(dN);