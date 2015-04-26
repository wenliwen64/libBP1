function skew=kurt2(x0)
% Getting the skewness (-90 to 90 degree) of profile x which a minimum
% Kurtosis
x0=x0-mean(x0);
n=181;
theta=linspace(-90,90,n)+k0;
kurtosis=zeros(1,n);
m=length(x0);

for i=1:n

        kurtosis(i)=mean(xx1.^4)/(mean(xx1.^2))^2-3; %kurtosis

end

[val index]=min(kurtosis);% finding the minimum;
skew=theta(index);

