%% I.Load data and do interpolation
close all;
clear all;
filename = '../logfile_2014_04_02_0_6';
logfile = importdata(filename);
load ../eqinfo_nchile_20100101_20140403_local.dat;
origin_power = logfile(:,4);
time_series = logfile(:,1);
ind_normal = find(origin_power ~= 0);
ind_missing = find(origin_power == 0);
sig_num = 1:numel(origin_power);
inter_power = interp1(ind_normal, origin_power(ind_normal), sig_num, 'nearest');

log10_inter_power = log10(inter_power);
%different from Sumatra
Mag = (log10_inter_power - 1.104) / 1.005;

%% II.Plot after moving average
figure(1);
hold all;
y = moving_average(Mag, 5);
plot(logfile(:,1), y, 'b-');
xlabel('time(s)');
ylabel('Magnitude');
title('comparison between interpolation(smooth) and ');

% Get the MAD spectrum
median_1 = median(y);
mad = median(abs(y-median_1));
mad_threshold = median_1 + mad;
%legend('mag\_original', 'mag\_after\_average');

%% III.Do peak finding
[PKS_y,LOCS_y] = findpeaks(y);
N = length(PKS_y);
kt = 1;
ts = 6;
tl = 40;
PKS0_y = zeros(numel(PKS_y), 1);
LOCS0_y = zeros(numel(PKS_y), 1);
stalta_ratio = zeros(numel(PKS_y), 1);
ratio_threshold = 1.165;

for i = 1:N
	 STA = median(y(max([1 LOCS_y(i)-ts/2]):min([length(logfile) LOCS_y(i)+ts/2])));
     LTA = median(y(max([1 LOCS_y(i)-tl/2]):min([length(logfile) LOCS_y(i)+tl/2])));
     %if PKS(i)>2*mean(Power0(max([1 LOCS(i)-40]):min([length(logfile) LOCS(i)+40])))
     stalta_ratio(i) = STA/LTA;
     if STA/LTA > ratio_threshold
        PKS0_y(kt) = PKS_y(i);
        LOCS0_y(kt) = LOCS_y(i);
        kt = kt+1;
     end
end

PKS0_y(PKS0_y == 0) = [];
LOCS0_y(LOCS0_y == 0) = [];
%catalogtime = datenum([eqinfo(:,1:6)]);
begin = datenum(2012,04,11,08,30,00);
time_bp_peak = logfile(LOCS0_y,1);
mag_bp_peak = y(LOCS0_y);
plot(time_bp_peak,mag_bp_peak,'r.','Markersize',20);
%TODO:plot((catalogtime-begin)*24*3600+520, eqinfo(:,9), 'g.', 'Markersize', 20);

%% IV.Old method to find out rising positions
figure(2);
hold on;
rising_position = zeros(numel(LOCS0_y),1);
inter_power_diff = diff(y);

% kt-1 is the number of peaks found using above method
for i = 1:(kt-1)
    for j = 1:LOCS0_y(i)
        if(inter_power_diff(LOCS0_y(i)-j) < 0.0)
           rising_position(i) = LOCS0_y(i) - j + 1;
           break;
        end
    end
end

mag_rising = y(rising_position);
time_rising = logfile(rising_position,1);
plot(time_bp_peak, mag_bp_peak, 'r.', 'Markersize', 25);
%TODO:plot((catalogtime-begin)*24*3600+520, eqinfo(:,9), 'g.', 'Markersize', 20);
plot(time_rising, mag_rising, 'k.', 'Markersize', 20);
plot(logfile(:,1), y, 'b-');
saveas(gcf, 'initiation_time_series', 'jpg');
print('-dpdf', '-r300', 'initiation_time_series.pdf');

%% V.Using trigger window to get the valley point position

tri_threshold = mad_threshold;    %3.312 This value is obtained by computing the MAD(median of absolute deviation)
diff_y = diff(y);
rising_y = zeros(1000,1);
count_rising_y = 0;
begin_win_array = zeros(1000,1);
end_win_array = zeros(1000,1);
count_win = 0;

for i = 1:numel(y)
     if count_win >= 1
         if (i <= end_win_array(count_win))
             continue;
         end
     end
    
    if y(i) >= tri_threshold
        begin_win = i;
        %search from current position to the end
        for j = i:numel(y)
            if y(j) <= tri_threshold
                end_win = j;
                break;
            end
        end
        
        %check if we can find real peaks within this triggered time window
        any_peak = 0;
        for k = begin_win:end_win
            if any(abs(k-LOCS0_y)==0)
                any_peak = 1;
                display('GOT ONE!');
                break;
            end
        end
        
        %save the qualified time window inoformation
        if any_peak == 1
            count_win = count_win + 1;
            begin_win_array(count_win) = begin_win;
            end_win_array(count_win) = end_win;
            display('SAVE ONE!');
        else
            continue; %continue i loop
        end
        
        %go back to search for the valley point
        if count_win == 1
            distance_from_last_win_end = begin_win_array(count_win) - 1; 
        else
            distance_from_last_win_end = begin_win_array(count_win) - end_win_array(count_win - 1); 
        end
        
        for k = 1:distance_from_last_win_end
            if diff_y(begin_win - k) < 0
                count_rising_y = count_rising_y + 1;
                rising_y(count_rising_y) = begin_win - k + 1;
                break;
            end
        end
        
    else
        continue;
    end
end 

begin_win_array(begin_win_array == 0) = [];
end_win_array(end_win_array == 0) = [];
%% VI.Plot the duration 
rising_y(rising_y == 0) = [];
figure(3);
hold all;
plot(logfile(:,1),y,'b-');
time_risingII = logfile(rising_y, 1);
mag_risingII = y(rising_y);
plot(time_risingII,mag_risingII,'k.','Markersize',20);
plot(time_bp_peak,mag_bp_peak,'r.','Markersize',20);
%plot((catalogtime-begin)*24*3600+520, eqinfo(:,9), 'g.', 'Markersize', 20);

%%=============TODO: we need to match the magnitude and then normalize the duration of the time and plot on the map
%By observing, we start from black dot 5 to 21.
count_dur = 0;
dur_pair = zeros(1000,2);

for i = 5:numel(rising_y)
    valley_p = rising_y(i);
                 
    %Find the first peak right after this valley point
    for j = 1:numel(LOCS0_y)
        diff_locs0 = LOCS0_y(j) - valley_p;
        if diff_locs0 >= 0
            peak_p = LOCS0_y(j);
            break;
        end
    end
    
    for j = valley_p:peak_p
        if y(j) < mad_threshold
            %display('continue');
            continue;
        else 
            %display('begin_end');
            count_dur = count_dur + 1;
            dur_pair(count_dur, 1) = j;
            dur_pair(count_dur, 2) = peak_p;
            y(j)
            break;
        end
    end
end

%dur_pair(dur_pair(:, 1) == 0) = [];

duration = zeros(count_dur, 1);
for i = 1:count_dur
    duration(i) = dur_pair(i, 2) - dur_pair(i, 1);
end


dur_corr_pair = [5, 10; 6, 12; 7, 13; 8, 17; 9, 21; 10, 22; 12, 25; 13, 26; 14, 27; 15, 28; 16, 30; 17, 31; 18, 34; 19, 35; 20, 36; 21, 37];
[num_pair, num_pair1] = size(dur_corr_pair);
matrix_af_norm = zeros(num_pair, 5);

for i = 1:count_dur
%     duration_id = dur_corr_pair(i, 1);
%     cat_id = dur_corr_pair(i, 2);
    peak_id_temp = dur_pair(i, 2);
    rising_id_temp = dur_pair(i, 1);
    mag_bp_temp = y(peak_id_temp);
    norm_factor = (mag_bp_temp / 6)^(1 / 3);   %(eqinfo(cat_id, 9)/6)^(1/3);
    lon_bp_temp = logfile(peak_id_temp, 3);
    lat_bp_temp = logfile(peak_id_temp, 2);
    duration_temp = duration(i);
    
    matrix_af_norm(i, 1) = lat_bp_temp;   %latitude
    matrix_af_norm(i, 2) = lon_bp_temp;   %logitude
    matrix_af_norm(i, 3) = mag_bp_temp;   %magnitude
    matrix_af_norm(i, 4) = duration_temp / norm_factor;%duration
    matrix_af_norm(i, 5) = duration_temp;
end

figure(4);
hold on;

scatter(matrix_af_norm(:,2), matrix_af_norm(:,1), 500*matrix_af_norm(:,3).^2/max(matrix_af_norm(:,3))^2, matrix_af_norm(:, 4), 'fill');
%M_usgs_total = horzcat(eqinfo(:, 8), eqinfo(:, 7), eqinfo(:,9)/max(eqinfo(:,9)),...
                       %(catalogtime-begin)*24);
%dlmwrite('total_usgs_sumatra', M_usgs_total, 'delimiter', ' ');
xlabel('lon(degree)');
ylabel('lat(degree)');
%plot(lon0, lat0, 'b*', 'Markersize', 20);
%plot(lon1, lat1, 'g*', 'Markersize', 20);
%title('usgs catalog');
box on;
daspect([1,cosd(2),1]);
colorbar;

%% VII.Check the duration vs. magnitude relationship
figure(5);
hold on;
peak_id = dur_pair(:, 2);
log10peak_dur = log10(y(peak_id(1:count_dur)));
plot(log(duration(1:count_dur)),  matrix_af_norm(:, 3), 'bo');