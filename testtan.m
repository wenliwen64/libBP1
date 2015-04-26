figure(1);
x=0:0.01:100;

for i=1:length(x)
c(i)=1*(1-atan(1*x(i))/(pi/2));
end
plot(x,c);