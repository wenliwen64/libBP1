

z=linspace(0,20,200);
for k=1:200
lhs1(k)=(1-1)*besselj(1,z(k))-z(k)*besselj(1+1,z(k));
lhs3(k)=(3-1)*besselj(3,z(k))-z(k)*besselj(3+1,z(k));
end
plot(z,lhs1,'r-',z,lhs3,'b-',z,0,'black-')
xlabel('Z');
ylabel('LHS');
legend('l=1','l=3');