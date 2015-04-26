function y=hyperbola(coef,x)
for k=1:1:length(x)
y(k)=-(coef(4)+coef(1)*x(k))/(coef(2)*x(k)+coef(3))
end
