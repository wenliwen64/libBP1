%% load data
clear
data=load('ARsample.txt');

y1=data(1:100,2);               % The first half data
%% AR prediction
figure;
plot(data(:,2),'k','linewidth',2); 
hold all;                       % plot the real data as a referance


for order=2:5:30                % let ar model order varies in 2:5:30 range
    
m=ar(y1,order,'burg');          % estimate ar model, y1*m = err white noise 

err=filter(m.a,m.c,y1);         % calculate err

err_predict=[err;err*0];        % Since we do not know future excitation error, 
                                % we just assumeit's zero
                                
y2(:,order)=filter(m.c,m.a,err_predict); % Make prediction

plot(y2(:,order));              % plot prediction

end

hold off

%%
% We can see the prediction series is very similar for difference AR orders
%% Spectrum
figure;
% real spectrum
amp=abs(fft(data(:,2)));
plot(amp(1:end/2),'k','linewidth',2); hold all
% padding zero spectrum
amp=abs(fft([y1;y1*0]));
plot(amp(1:end/2),'b','linewidth',2);
% AR spectrum
for order=[2,12,22]
amp=abs(fft([y2(:,order)]));
semilogy(amp(1:end/2));
end
hold off
%%
% We can see AR spectrum has more sharp peak than padding zeros spectrum
% and real spectrum, and also has small size ripple. So by extracting extra
% information from data, AR method can increase frequency resolution in
% some way.