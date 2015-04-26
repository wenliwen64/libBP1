clear;
a=6.6;
b=3.7;
z=0;
%c=roots([1/b^6 0 -8/b^4 0 (24/b^2-16/a^2) 0 -16*(1-b^2/a^2)])

c=3.4143;
sa=(1-c^2/a^2)^0.5;
sb=(1-c^2/b^2)^0.5;
t=linspace(0,9,201)
%c=roots([1/b^6 0 -8/b^4 0 (24/b^2-16/a^2) 0 -16*(1-b^2/a^2)])
f=1/b^6*t.^6-8/b^4*t.^4+(24/b^2-16/a^2)*t.^2-16*(1-b^2/a^2)
plot(t,f)
%u=(exp(1/c*sa*z)-(1-0.5*c^2/b^2)*exp(1/c*sb*z))*exp(i*t);
%w=-i*(sa*exp(1/c*sa*z)+1/sb*(1-0.5*c^2/b^2)*exp(1/c*sb*z))*exp(i*t);
%hold on;
%plot(u,w,'r')
xlabel('c');
ylabel('LHS');
%legend('z=0','z=1','z=3','z=4','z=5','z=8')
%title('plot of u vs w at w=1 and x=0')
