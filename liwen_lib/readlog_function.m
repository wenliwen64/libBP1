function readlog_function(file_name)
close all;
%%%=====================================%%%
logfile = importdata(file_name);
dates = sscanf(file_name, 'logfile_%d_%d_%d_%d_%d');
file_id =[num2str(dates(1)), '_', num2str(dates(2)), '_', num2str(dates(3)), '_', num2str(dates(4)), '_', num2str(dates(5))];
%logfile = importdata('./0831logfile_0_7.txt');
%load eqinfo.dat;
eqinfo = importdata('eqinfo_2014.dat');
lon0 = -70.817;% used for bp
lat0 = -19.642;% used for bp
% lon1=92.452;
% lat1=0.733;
dlon=-0.091; % catalog - dlon = lon0
dlat=0.07; % catalog - dlat = lat0
%%%=====================================%%%

%%=======================Figure1:comparison between catalogs========================

figure(1);     

%%%=============================%%%

Mag = log10(logfile(:, 4)); %1.4*((8.2-5.7)/(7.92-5.81))

%%%=============================%%%

Power0 = logfile(:, 4);
% Power0 = padarray(Power0, [0 20], 'pre');
% Power0 = padarray(Power0, [0 19], 'post');
Power = log10(logfile(:, 4));
% Power = padarray(Power, [0 20], 'pre');
% Power = padarray(Power, [0 19], 'post');

[PKS, LOCS] = findpeaks(Power0);
N = length(PKS);
kt = 1;
ts = 4;
tl = 60;
for i = 1:N
	 STA = median(Power0(max([1 LOCS(i)-ts/2]):min([length(logfile) LOCS(i)+ts/2])));
     LTA = median(Power0(max([1 LOCS(i)-tl/2]):min([length(logfile) LOCS(i)+tl/2])));

if STA/LTA > 2.7
        PKS0(kt)=PKS(i);
        LOCS0(kt)=LOCS(i);
        kt = kt + 1;
     end
end
figure(1);
hold on;

%%%===================================================================%%%
catalogtime = datenum([eqinfo(:,1:6)]);
%%%===================================================================%%%

plot(logfile(LOCS0,3),logfile(LOCS0,2),'r.','Markersize',20);
%begin=datenum(2012, 04, 11, 08, 30,00);
%over=datenum([2014 04 02 5 0]);
%I=find(eqinfo(:,7)>begin&eqinfo(:,7)<over);

%%%====================Set up the time shift to align the time series==============================%%%

time_span = 2.0/24;

this_file_begin_utc = datenum([dates(1) dates(2) dates(3) dates(4) 0 0]);
traveltime = 0.008229167; %datenum of travel time calculated by the Mw 8.0 earthquake
I = find(eqinfo(:, 7) >= this_file_begin_utc - traveltime & eqinfo(:,7) < this_file_begin_utc - traveltime + time_span);
plot(eqinfo(I,9) - dlon, eqinfo(I,8) - dlat, 'b.', 'Markersize', 20);
legend('BP aftershock', 'Chile Local catalog');
II = find(eqinfo(:,10) > 8.0); %find out the main shock
%%%==================================================%%%

	
%plot(lon0,lat0,'b*','Markersize',20);
%plot(lon1,lat1,'g*','Markersize',20);
ylabel('Latitude (^o)');
xlabel('Longitude (^o)');
%box on;
%daspect([1, cosd(2),1]);
% 
% for i = 1:numel(I)
%     text(eqinfo(LOCS0(i), 3), logfile(LOCS0(i), 2), num2str(i));
% end
print('-dpdf','-r300','BPaftershock_nchile.pdf'); 

saveas(gcf, ['BPaftershock_nchile', '_', file_id], 'jpg');

%%=======================Figure2: backprojection time series===============
figure(2);
hold on;
plot(logfile(LOCS0,1)/3600, Mag(LOCS0), 'r.', 'Markersize',20);
% for i = 1:numel(LOCS0)
%    line([logfile(LOCS0(i),1) / 3600 logfile(LOCS0(i),1) / 3600], [0 Mag(LOCS0(i))]);
% end

%%%====================================================================%%%
plot((eqinfo(I,7)-this_file_begin_utc+traveltime)*24, eqinfo(I,10),'g*','Markersize',10);
%%%====================================================================%%%

% for i = 1:numel(I)
%    line([(eqinfo(I(i),7) - begin)*24+0.198 (eqinfo(I(i),7) - begin)*24+0.198], [0 eqinfo(I(i),11)]);
% end
xlabel('time (sec)');
ylabel('Magnitude');
legend('BP catalog','USGS catalog');
plot(logfile(:, 1)/3600, Mag,'-');
box on;
%daspect([1, cosd(2),1]);
for i = 1:numel(LOCS0)
    text(logfile(LOCS0(i), 1)/3600, Mag(LOCS0(i)), num2str(i));
end

for i = 1:numel(I)
    text((eqinfo(I(i), 7)-this_file_begin_utc+traveltime)*24, eqinfo(I(i), 10), num2str(i));
end


print('-dpdf','-r300','BPcatalogTimeSeries_nchile.pdf'); 
saveas(gcf, ['BPcatalogTimeSeries_nchile', file_id], 'jpg');



%============points in catalog having corresponding pb ponits, picked by hand===========
ind_cat = zeros(1, 100);
ind_bp = zeros(1, 100);

cor_count = 0;
for i=1:numel(I)
    time_eq = (eqinfo(I(i), 7)-(this_file_begin_utc-traveltime))*24;  % in hours
    if(min(abs(time_eq - logfile(LOCS0, 1)/3600)) > 100/3600)
        continue;
    else
        cor_count=cor_count+1;
        [Y, loc_min] = min(abs(time_eq - logfile(LOCS0, 1)/3600));
        ind_cat(cor_count)=i;
        ind_bp(cor_count)=loc_min;
        
    end
end

ind_cat(ind_cat==0)=[];
ind_bp(ind_bp==0)=[];
% ind_cat = [2 6 7 9 10 12 16 20 21 28];
% ind_bp =  [2 3 4 7 11 15 16 17 19 22];
%cor_cat = [1 2 3 4 6 7 8 9 10 11 13 14 16 17 18 20 22 24 25 26 27 30  31 34];
% III = zeros(numel(cor_cat), 1);
% Y = zeros(numel(LOCS0), numel(cor_cat));
% 
% for i = 1:numel(cor_cat)
%     [KKKK, III(i, 1)] = min(abs(logfile(LOCS0,1)/3600 - ((eqinfo(I(i),7) - begin)*24+0.198)));
%     Y(:,i) = abs(logfile(LOCS0,1)/3600 - ((eqinfo(I(i),7) - begin)*24+0.198));
% 
% end
% III

%%======================Figure3: back-projection spatial-temporal map===================
% figure(3);
% hold on;
% %%%=============================================%%%
% plot(eqinfo(II, 9), eqinfo(II, 8), 'b*', 'Markersize', 40);
% %%%=============================================%%%
% %plot(lon1,lat1,'g*','Markersize',60);
% scatter(logfile(LOCS0,3), logfile(LOCS0,2), 100*log10(logfile(LOCS0,4))/max(log10(logfile(LOCS0,4))), (LOCS0-LOCS0(1))/3600,'fill');
% xlabel('lon(degree)');
% ylabel('lat(degree)');
% title('bp result');
% colorbar;
% box on;
% daspect([1, cosd(2), 1]);
% print('-dpdf','-r300','BP_SpatialTemporal.pdf'); 
% saveas(gcf, 'BP_SpatialTemporal', 'jpg');
% 
% %%===================Figure4: usgs patial-temporal map======================
% figure(4);
% hold on;
% 
% %%%==========================================================================%%%
% catalog_time_concerned = eqinfo(I, 7);
% scatter(eqinfo(I, 9),eqinfo(I, 8), 100*eqinfo(I, 11)/max(eqinfo(I,11)),(catalog_time_concerned - begin) * 24,'fill');
% %%%===========================================================================%%%
% 
% xlabel('lon(degree)');
% ylabel('lat(degree)');
% %plot(lon0,lat0,'b*','Markersize',20);
% 
% %%%==================================================================%%%
% plot(eqinfo(II, 9), eqinfo(II, 8), 'b*', 'Markersize', 40);
% title('usgs catalog');
% %%%===============================================================%%%
% 
% box on;
% daspect([1,cosd(2),1]);
% colorbar;
% print('-dpdf','-r300','ChileLocalCatalog_SpatialTemporal.pdf'); 
% saveas(gcf, 'ChileLocalCatalog_SpatialTemporal', 'jpg');

%%====================Figure5: connect back-projection and usgs corresponding aftershocks===================================================
figure(5);
hold on;
% ind_bp = III;
%ind_usgs = [1 2 3 4 6 7 9 10 12 13 17 21 22 24 25 27 28 30 31 35 36 37 39];
ind_usgs = ind_cat;

lat_usgs = eqinfo(I(ind_usgs), 8) - dlat;
lon_usgs = eqinfo(I(ind_usgs), 9) - dlon;

lat_bp = logfile(LOCS0(ind_bp), 2);
lon_bp = logfile(LOCS0(ind_bp), 3);

lat_diff = lat_bp - lat_usgs;
lon_diff = lon_bp - lon_usgs;
quiver(lon_usgs, lat_usgs, lon_diff, lat_diff);
%quiver(zeros(size(lon_diff)), zeros(size(lon_diff)), lon_diff, lat_diff);
catalogtime = eqinfo(ind_usgs,7);
%scatter(eqinfo(:,8),eqinfo(:,7), 100*eqinfo(:,9)/max(eqinfo(:,9)),(catalogtime-begin)*24,'fill');
%scatter(eqinfo(ind_usgs, 8), eqinfo(ind_usgs,7), 100*eqinfo(ind_usgs,9)/max(eqinfo(ind_usgs,9)), 'b','fill');%(catalogtime-begin)*24,'fill');
%colorbar;
%scatter(logfile(ind_bp,3), logfile(ind_bp,2), 100*log10(logfile(ind_bp,4))/max(log10(logfile(ind_bp,4))), 'r', 'fill');%(LOCS0(ind_bp))/3600,'fill');
%colorbar;
%plot(logfile(LOCS0,3),logfile(LOCS0,2), 'r.','Markersize',20);
plot(lon_bp, lat_bp, 'r.','Markersize',20);
plot(lon_usgs, lat_usgs, 'b.','Markersize',20);

%plot(eqinfo(I,9),eqinfo(I,8), 'b.','Markersize',20);
legend('shift','bp','usgs');
xlabel('lon(degree)');
ylabel('lat(degree)');
%plot(lon0,lat0,'b*','Markersize', 20);
%plot(lon1,lat1,'g*','Markersize', 20);
for i = 1:numel(ind_bp)
    plot([lon_bp(i) lon_usgs(i)], [lat_bp(i) lat_usgs(i)], 'k-');
end
daspect([1,cosd(2),1]);
box on;
print('-dpdf','-r300','BP_CHILELOCAL_Shift.pdf');
saveas(gcf, ['BP_CHILELOCAL_Shift', file_id], 'jpg');

%%=========Figure 6:bp result to get the path of propagation =============
% figure(6);
% hold on;
% plot(lon0,lat0,'b*','Markersize',60);
% plot(lon1,lat1,'g*','Markersize',60);
% scatter(logfile(LOCS0,3), logfile(LOCS0,2), 100*log10(logfile(LOCS0,4))/max(log10(logfile(LOCS0,4))), (LOCS0-LOCS0(1))/3600, 'fill');
% xlabel('lon(degree)');
% ylabel('lat(degree)');
% title('bp result');
% colorbar;
% box on;
% daspect([1,cosd(2),1]);
% 
% [X, Y] = ginput(2);
% plot(X(1), Y(1), 'r^', 'Markersize', 10);
% plot(X(2), Y(2), 'r^', 'Markersize', 10);
% line(X,Y);
% rect_x = 90.8;
% rect_y = 0.1;
% rect_w = 1.8;
% rect_h = 1.8;
% display('straight line segment position');
% X
% Y
% text(X(1), Y(1), ['(' num2str(X(1)) ',' num2str(Y(1)) ')']);
% text(X(2), Y(2), ['(' num2str(X(2)) ',' num2str(Y(2)) ')']);
% rectangle('position', [rect_x, rect_y, rect_w, rect_h]);
% print('-dpdf','-r300','BP_SpatialTemporal_path.pdf');
% saveas(gcf, 'BP_Spatial_Tmeporal_path', 'jpg');

%%===========Figure 7: rupture path=================================
% figure(7);
% slope = (Y(2) - Y(1))/(X(2) - X(1));
% point1 = [X(1), Y(1)];
% point2 = [X(2), Y(2)];
% 
% event_count = 1;
% event_location = zeros(1000, 2);
% event_time = zeros(1000, 1);
% event_mag = zeros(1000,1);
% for i = 1:numel(LOCS0)
%     %logfile(LOCS0,3), logfile(LOCS0,2);
%     if(logfile(LOCS0(i), 3) < rect_x | logfile(LOCS0(i), 3) > rect_x + rect_w) 
%         continue;
%     end
%     if(logfile(LOCS0(i), 2) < rect_y | logfile(LOCS0(i), 2) > rect_y + rect_h) 
%         continue;
%     end
%     event_location(event_count, 1) = logfile(LOCS0(i), 3); % first column is logitude
%     event_location(event_count, 2) = logfile(LOCS0(i), 2); % second column is latitude
%     event_time(event_count, 1) = logfile(LOCS0(i), 1);
%     event_mag(event_count, 1) = logfile(LOCS0(i), 4);
%     event_count = event_count + 1; 
% end
% 
% distance_from_point1 = zeros(1, event_count-1);
% distance_from_line = zeros(1, event_count-1);
% 
% 
% for i = 1:event_count-1
%     distance_from_point1(i) = abs(dot([event_location(i,1),event_location(i,2)] - point1, point2 - point1))/sqrt(dot(point2 - point1, point2 - point1));
%     
%     length_from_point1 = dot([event_location(i,1),event_location(i,2)] - point1, [event_location(i,1),event_location(i,2)] - point1);
%     distance_from_line(i) = sqrt(length_from_point1 - distance_from_point1(i)*distance_from_point1(i));
% end
% scatter(event_time(1:(event_count-1)), distance_from_point1/180*pi*6373, 100*log10(event_mag(1:(event_count-1)))/max(log10(event_mag(1:(event_count)))), distance_from_line, 'fill');
% colorbar;
% distance_from_line
% title('sumatra Apr 2014 rupture distance vs. time');
% xlabel('time(second)');
% ylabel('rupture path distance(km)');

%%===============Figure 8. rising points in time series=========================
% figure(8);
% hold on;
% rising_position = zeros(numel(LOCS0),1);
% amp_array = logfile(:, 4);
% amp_array_der = diff(amp_array);
% for i = 1:numel(LOCS0)
%     for j = 1:LOCS0(i)
%         if(amp_array_der(LOCS0(i)-j) < 0.0)
%            rising_position(i) = LOCS0(i) - j + 1;
%            
%            break;
%         end
%     end
% end
% catalogtime = datenum([eqinfo(:,1:6)]);
% 
% M_rising_points = horzcat(logfile(LOCS0, 1), logfile(LOCS0, 2), logfile(LOCS0, 3), ...
%                          logfile(rising_position, 1), logfile(rising_position, 2), logfile(rising_position, 3), logfile(LOCS0, 4));
% dlmwrite(['rising_points_data_', file_id, '.txt'], M_rising_points, 'delimiter', ' ');
figure(8)
hold on;
rising_position = zeros(numel(LOCS0),1);
amp_array = logfile(:, 4);
amp_array_der = diff(amp_array);
LOCS0_temp = LOCS0;

for i = 1:numel(LOCS0)
    for j = 1:LOCS0(i)
        if(amp_array_der(LOCS0(i)-j) < 0.0)
           rising_position(i) = LOCS0(i) - j + 1;
           
           break;
        end
    end
end
mag_m = logfile(:, 4);
mag_m(mag_m<=0)=[];
mean_mag = mean(log10(mag_m));
for i = 1:numel(rising_position)
    peak_mag = log10(logfile(LOCS0(i), 4));
    valley_mag = log10(logfile(rising_position(i), 4));
    peak_mag - mean_mag
    abs(valley_mag - mean_mag)
    if(abs(valley_mag - mean_mag) > (peak_mag - mean_mag)/ 2.0)
        rising_position(i) = 0;
        LOCS0_temp(i) = 0;
    end 
end
mean_mag
rising_position
LOCS0_temp
rising_position(rising_position == 0) = [];
LOCS0_temp(LOCS0_temp == 0) = [];
rising_position
LOCS0_temp
catalogtime = datenum([eqinfo(:,1:6)]);

M_rising_points = horzcat(logfile(LOCS0_temp, 1), logfile(LOCS0_temp, 2), logfile(LOCS0_temp, 3), ...
                         logfile(rising_position, 1), logfile(rising_position, 2), logfile(rising_position, 3), log10(logfile(LOCS0_temp, 4)));
dlmwrite(['rising_points_data_', file_id, '.txt'], M_rising_points, 'delimiter', ' ');
% plot(logfile(LOCS0,1),Mag(LOCS0),'r.','Markersize',20);
% plot((catalogtime-begin)*24*3600+520,eqinfo(:,9),'g*','Markersize',20);
% plot(logfile(rising_position,1),Mag(rising_position),'k.','Markersize',20);
% plot(logfile(:,1),Mag,'-');
% saveas(gcf, 'initiation_time_series', 'jpg');
% print('-dpdf','-r300','initiation_time_series.pdf');

%%================Figure 9. rising points on map==========================
% figure(9);
% hold on;
% scatter(logfile(LOCS0,3),logfile(LOCS0,2), 100*log10(logfile(LOCS0,4))/max(log10(logfile(LOCS0,4))),'fill');
% plot(logfile(rising_position,3),logfile(rising_position,2),'r.','Markersize',20);
% 
% lat_bp = logfile(LOCS0, 2);
% lon_bp = logfile(LOCS0, 3);
% 
% lat_rising = logfile(rising_position, 2);
% lon_rising = logfile(rising_position, 3);
% 
% lat_diff = lat_rising - lat_bp;
% lon_diff = lon_rising - lon_bp;
% quiver(lon_bp, lat_bp, lon_diff, lat_diff);
% plot([lon_bp lon_rising]', [lat_bp lat_rising]', 'k-');
% saveas(gcf, 'initiation_vs_peak_map', 'jpg');
% print('-dpdf','-r300','initiation_vs_peak_map.pdf');

%%===============Figure 10. bp vs. catalog with magnitude===============
figure(10);
hold on;
scatter(logfile(LOCS0,3), logfile(LOCS0,2), 500*log10(logfile(LOCS0,4)).^2/max(log10(logfile(LOCS0,4)))^2, 'r', 'fill');
scatter(eqinfo(I,9) - dlon, eqinfo(I,8) - dlat, 500*eqinfo(I,10).^2/max(eqinfo(I,10))^2,'b');
plot(eqinfo(II,9), eqinfo(II,8), 'g*', 'Markersize', 20);
%quiver(lon_usgs, lat_usgs, lon_diff, lat_diff);
for i = 1:numel(ind_bp)
    plot([lon_bp(i) lon_usgs(i)], [lat_bp(i) lat_usgs(i)], 'k-');
end
legend('bp', 'chile\_catalog');
title('bp vs. chile\_catalog');
saveas(gcf, ['bp_vs_catalog_mag_space', file_id], 'jpg');

%%===============Figure 11. ===================================
shift_bp_usgs = sqrt(lat_diff .* lat_diff + lon_diff .* lon_diff);
az_bp_usgs = mod(atan2(lat_diff, lon_diff), 2*pi);
M_bp_usgs = horzcat(shift_bp_usgs, az_bp_usgs, log10(logfile(LOCS0(ind_bp), 4)), eqinfo(I(ind_usgs), 10), logfile(LOCS0(ind_bp), 3), logfile(LOCS0(ind_bp), 2) , eqinfo(I(ind_usgs), 9)-dlon, eqinfo(I(ind_usgs), 8)-dlat);
dlmwrite(['shift_bp_usgs_', file_id, '.txt'], M_bp_usgs, 'delimiter', ' ');



%%==============Figure 12. Magnitude Relationship================
figure(11);
scatter(Mag(LOCS0(ind_bp)), eqinfo(I(ind_cat), 10));

%%==============Record information===============================
M_bp_sum = horzcat(logfile(LOCS0, 3), logfile(LOCS0, 2), log10(logfile(LOCS0, 4)));
dlmwrite(['bp_peak_data_', file_id, '.txt'], M_bp_sum, 'delimiter', ' ');
M_cat_sum = horzcat(eqinfo(I, 8) - dlat, eqinfo(I, 9) - dlon, eqinfo(I, 10));
dlmwrite(['cat_peak_data_', file_id, '.txt'], M_cat_sum, 'delimiter', ' ');
end