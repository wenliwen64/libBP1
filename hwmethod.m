clear;
load i_p_delta_2;
load jeffreys_model;
dist=i_p_delta_2(:,1);
t=i_p_delta_2(:,2);
p=i_p_delta_2(:,3);
r0=6371;
figure;
for j=2:length(dist)
    pc(j)=(t(j)-t(j-1))/(dist(j)-dist(j-1))/pi*180;
end
plot(dist,t/60,'r-',dist,p/60,'b-',dist(2:50),pc(2:50)/60,'g-');
title('plot of t and p as a function of distance');
xlabel('distance (deg)');
ylabel('t and p (min)');
legend('trave time','ray parameter','p caculated by differentiating t');
for k=1:length(dist)
    integral=0;
    for kk=2:k
    integral=integral+(log(p(kk)/p(k)+((p(kk)/p(k))^2-1)^0.5)+log(p(kk-1)/p(k)+((p(kk-1)/p(k))^2-1)^0.5))/180*pi;
    end
    r1(k)=r0*exp(-1/pi*integral);
    v(k)=r1(k)/p(k);
    depth(k)=r0-r1(k);
end


figure;
plot(depth,v,'b',jeffreys_model(:,1),jeffreys_model(:,2),'r');
title('plot of velocity as a function depth');
xlabel('depth (km)');
ylabel('velocity (km)');
legend('Herglotz-Wiechert Method','Jeffreys model');

