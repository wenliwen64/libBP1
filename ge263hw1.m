h=1;
D=1;
dt=1*h^2/D;
t=1:100;
wx=t*0;
for j=45:1:50
    wx(j)=1;
end
for j=1:100
    for k=1:100
        if j==k
            T(j,k)=-2;
        else if abs(j-k)==1
                T(j,k)=1;
            else
                T(j,k)=0;
            end
        end
    end
end
for j=1:100
    for k=1:100
        if j==k
            I(j,k)=1;
       
        else
            I(j,k)=0;
        end
            
    end
end

figure(1);
%title('Heat flow example for Finite difference dt=1/4*h^2/D')

subplot(4,2,1);
Q=wx;
for j=1:1
    %Q=Q*(I+D*dt/h^2*T);
    %Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('explicit 0');
ylim([0 1.1]);

subplot(4,2,2);
Q=wx;
for j=1:100
    Q=Q*(I+D*dt/h^2*T);
    %Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('explicit 100');
%ylim([0 1.1]);


subplot(4,2,3);
Q=wx;
for j=1:1000
    Q=Q*(I+D*dt/h^2*T);
    %Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('explicit 1000');
%ylim([0 1.1]);


subplot(4,2,4);
Q=wx;
for j=1:2000
    Q=Q*(I+D*dt/h^2*T);
    %Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('explicit 2000');
%ylim([0 1.1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%title('Heat flow example for Finite difference dt=1/4*h^2/D')

subplot(4,2,5);
Q=wx;
for j=1:1
    %Q=Q*(I+D*dt/h^2*T);
    %Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('implicit 0');
ylim([0 1.1]);

subplot(4,2,6);
Q=wx;
for j=1:100
    %Q=Q*(I+D*dt/h^2*T);
    Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('implicit 100');
%ylim([0 1.1]);


subplot(4,2,7);
Q=wx;
for j=1:1000
    %Q=Q*(I+D*dt/h^2*T);
    Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('implicit 1000');
%ylim([0 1.1]);


subplot(4,2,8);
Q=wx;
for j=1:2000
    %Q=Q*(I+D*dt/h^2*T);
    Q=Q*(I+D*dt/h^2*T)*inv(I-D*dt/h^2*T);
end
plot(Q);
xlabel('x');
ylabel('Q');
title('implicit 2000');
%ylim([0 1.1]);



