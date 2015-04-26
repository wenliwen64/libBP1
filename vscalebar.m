function vscalebar(x0,y,Y1,Y2,C,shift,linewidth)
%plot a vertical scale bar, 
%x0: x location of the bar
%y: y location of the top and bottom of the bar
%Y1, Y2: number to be plotted on top and bottom pf the bar
%C: color
%shift: units in x to shift the number
%linewidth: linewidth

y1=y(1);
y2=y(2);
shift1=shift(1);
shift2=shift(2);
w=0.1*(y2-y1);
plot([x0 x0],[y1 y2],C,'LineWidth',linewidth);
plot([x0-w x0+w],[y1 y1],C,'LineWidth',linewidth);
plot([x0-w x0+w],[y2 y2],C,'LineWidth',linewidth);
text(x0-shift1,y1,Y1,'color',C);
text(x0-shift2,y2,Y2,'color',C);