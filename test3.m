close all;

hold on;
xx=linspace(0,30,6);
ddx=ones(1,6);
       for dd=10:5:30
           
           
           plot((xx.^2+dd^2).^0.5,ddx*dd,'^-');
           for cc=1:length(xx)
           text((xx(cc)^2+dd^2).^0.5,dd+1,num2str(xx(cc)));
           end
       end