clear;
AU=1.49598e11
GM=1.9891e30*6.67e-11
a =linspace(0,-1.1,12)*1/AU/1000
v0=30e3
for k=1:11
    NN(k)=-(a(k+1)-a(k))*GM/(2*v0)*exp(1/v0^2*a(k+1)*GM)
end
bar(a(1:11)*AU*1000,NN*26/sum(NN));
ylim([2.36 2.365])
    

