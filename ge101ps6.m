clear;
x=-90:10:90;
y=htr(x);


l=1e15*[0 -0.4 -1.5 -3.1 -4.6 -5.5 -5.3 -4.3 -2.4 0.3 2.9 4.9 5.6 5.4 4.6 3.2 1.7 0.5 0];
plot(x,y,'r',x,l,'b');
xlabel('Lattitude(from South Pole to North Pole)');
ylabel('Net Transport Heat(W)')
Title('Poleward transport of energy(red for model, blue for observed data)')
figure;
%t(2)=1;
for k=1:length(l)-1
    df(k)=(l(k+1)-l(k))/10
end
    

for k=1:length(l)-1
    t(k)=tef(df(k),x(k)+5)
end

plot(x(1:18)+5,t,'b');
xlabel('Lattitude(from South Pole to North Pole)');
ylabel('effective radiating temperature Te')
Title('Te as a function of latitude')