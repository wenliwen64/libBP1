function ydata=finddipping(x,xdata)
N=length(xdata)/2;
str=x(1);
dip=x(2);
v1=2.0;
v2=1;
for i=1:N
    tbaz(i)=xdata(i);
    slowt(i)=xdata(i+N);
    [daz,dsl]=azsldipping(str,dip,tbaz(i),slowt(i),v1,v2);
    ddaz(i)=daz;
    ddsl(i)=dsl;
end
ydata=[ddaz 100*ddsl];