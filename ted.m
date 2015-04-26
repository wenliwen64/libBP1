clear;
age=[9320 9200 8800 8000 6000 0]
w=[446 437 385 331 261 180]
C=10;
n=1;
for k=2 :6
    aaa(k-1)=((w(k)-w(k-1))/(age(k)-age(k-1))/C)^(1/n);
end
aaa
plot(aaa,'*-')