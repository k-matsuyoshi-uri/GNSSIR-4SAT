% compare GNSS-IR RH with RH from Tide Gauge MSL for 2 Constellations
% Kayli Matsuyoshi

%% setup 
close all;
clear;

%% user-defined parameters 

% see station_Metadata.m to create parameters
load('N300_Station_Metadata.mat'); % change to appropriate metadata file

% daily reflector height (RH) statistics directory saved from LSP analysis
SNR_stats_directory = 'SNR_stats/N300/'; 
    % files should be in SSSS_YYYY_MM_DD_SNR_stats.mat format
    
% reflector height (RH) txt data directory saved from LSP analysis
RH_txt_directory = 'RH_S1/N300/';

% number of constellations observed
num_constellations = 2;

%% tide gauge adjustment distance

% distance between tidal datum epoch MSL and GNSS antenna (m)
    % Ben's Method: 5.0952 m (APC_to_MWWL + MWWL_to_STND - MSL_to_STND)
    % NOAA's Method: (APC_RTP + ARP_to_STND - MSL_to_STND)
GNSS_MSL_Distance = ARP_to_RLR - STND_to_RLR - MSL_to_STND + APC_to_ARP;

%% load RH daily statistics

% SNR_stats/ files contain SNR data, SNR statistics, and RH statistics
    % here only the RH statistics are extracted
snr_meta = struct2table(dir(SNR_stats_directory)); 
files = (snr_meta(3:end,1)); % first two rows are folder data, not a file
median_RH = nan(size(files,1), 1);
mean_RH = median_RH; stdev_RH = median_RH; day_GNSS = median_RH;

for f = 1:size(files,1)
    filename =  string(table2cell(files(f,1)));
    load(append(SNR_stats_directory, filename)); 
    % files should be in SSSS_YYYY_MMDD_SNR_stats.mat format
        % SSSS = station name
        % YYYY = year
        year = str2double(eraseBetween(eraseBetween(filename, 1, 5), 5, 23));
        % MM = month
        month = str2double(eraseBetween(eraseBetween(filename, 1, 10), 3, 18));
        % DD = day
        day = str2double(eraseBetween(eraseBetween(filename, 1, 12), 3, 16));
    % save avg RH, median RH, and stdev of RH for each day
    day_GNSS(f) = datenum(year, month, day);
    median_RH(f) = RH.med_RH;
    mean_RH(f) = RH.avg_RH;
    stdev_RH(f) = RH.std_dev_RH;
end

%% load and combine RH measurements

% RH_S1/ files contain LSP outputs of RH for each satellite that passes
% quality control
    % once the txt file data is concatenated, a time series can be made
rh_meta = struct2table(dir(RH_txt_directory));
files = (rh_meta(3:end,1)); % first two rows are folder data, not a file
rh_table = []; % expand matrix of data when concatenating
for f = 1:size(files,1)
    filename = string(table2cell(files(f,1)));
    rh_table = [ rh_table; load(append(RH_txt_directory, filename))]; %#ok<*AGROW>
end

% extract time data
    % column 13 contains UTC hour in decimal
time_GNSS = datenum(rh_table(:,1), rh_table(:,2), rh_table(:,3)) + (rh_table(:,13)/24); 

% make the entire matrix chronological
[time_GNSS, chronological] = sort(time_GNSS);
rh_table = rh_table(chronological, :);
RH_GNSS = rh_table(:,4);

% index by satellite constellation
sat_num = rh_table(:,7);
gps = find(sat_num < 100);
glonass = find(100 <= sat_num & sat_num < 200);

stdev_RH_timeseries = nan(size(time_GNSS));
% assign daily stdev_RH for each time stamp for RMSE calculation
for i = 1:length(day_GNSS)
    if i < length(day_GNSS)
        stdev_RH_timeseries(time_GNSS >= day_GNSS(i) & time_GNSS < day_GNSS(i+1)) = stdev_RH(i);
    else
        stdev_RH_timeseries(time_GNSS >= day_GNSS(i)) = stdev_RH(i);
    end
end

%% load TG MSL data

% TG data contains hourly MSL from NOAA Tides and Currents
warning off;
TG = readtable(TG_MSL_data); warning on;
time_TG = datenum(append(TG.Date,' ', TG.Time_GMT_));
% calculate a mean sea level for the whole dataset (m)
MSL = mean(TG.Verified_m_, 'omitnan');

%% daily mean TG

% collect daily average TGA to compare to daily GNSS-IR RH 
TG_daily_mean = nan(size(day_GNSS));
for i = 1:length(day_GNSS)
   if i < length(day_GNSS)
       TG_daily_mean(i) = mean(TG.Verified_m_(time_TG >= day_GNSS(i) & time_TG < day_GNSS(i+1)));
   else
       TG_daily_mean(i) = mean(TG.Verified_m_(time_TG >= day_GNSS(i)));
   end
end

%% compute sea level from GNSS-IR measurements

TG_daily_mean = TG_daily_mean(~isnan(TG_daily_mean));
day_GNSS = day_GNSS(~isnan(TG_daily_mean));
stdev_RH = stdev_RH(~isnan(TG_daily_mean));

% daily median mean sea level
SL_GNSSIR_MED = GNSS_MSL_Distance - median_RH;
SL_GNSSIR_MED = SL_GNSSIR_MED(~isnan(TG_daily_mean));

% daily average mean sea level
SL_GNSSIR_MEAN = GNSS_MSL_Distance - mean_RH;
SL_GNSSIR_MEAN = SL_GNSSIR_MEAN(~isnan(TG_daily_mean));

% conversion of GNSS-IR RH time series to GNSS-IR sea level time series
SL_GNSSIR = GNSS_MSL_Distance - RH_GNSS;

% interpolate tide gauge sea level time series to GNSS-IR time series
SL_TG = interp1(time_TG, TG.Verified_m_, time_GNSS);

% remove empty values from unequal start and end times
SL_TG = SL_TG(~isnan(SL_TG));
SL_GNSSIR = SL_GNSSIR(~isnan(SL_TG));
SL_stdev = stdev_RH_timeseries(~isnan(SL_TG));
time_SL = time_GNSS(~isnan(SL_TG)); 
SL_TG = SL_TG(abs(SL_GNSSIR)<1);
SL_stdev = SL_stdev(abs(SL_GNSSIR)<1);
time_SL = time_SL(abs(SL_GNSSIR)<1);
SL_GNSSIR = SL_GNSSIR(abs(SL_GNSSIR)<1);

sat_index = struct;

% interpolated indices
sat_index.gps = gps(gps<=length(time_SL));
sat_index.glonass = glonass(glonass<=length(time_SL));

% combination indices 
sat_index.gps_glo = union(sat_index.gps, sat_index.glonass);

%% evaluations
% re-run section if needed

[evals, evaluations] = compute_Evaluations(num_constellations, SL_GNSSIR_MED, SL_GNSSIR_MEAN, SL_stdev, TG_daily_mean, SL_GNSSIR, SL_TG, sat_index);

% list of evaluations
constellations = {'GPS'; 'GLONASS'};
combinations_2 = {'GPS+GLO'};
extra = {'All Satellites'; 'Daily Average'};

all_labels = [constellations; combinations_2; extra];

% save evaluations and labels
%save([GNSS_station '_evaluations.mat'], 'evaluations', 'all_labels', 'evals', 'sat_index', 'SL_GNSSIR', 'SL_TG', 'SL_GNSSIR_MED', 'SL_GNSSIR_MEAN', 'TG_daily_mean', 'SL_stdev', 'time_SL');

%% temporal resolution calculation

% time difference between all GNSS-IR measurements in minutes
time_diff = (time_GNSS(2:end) - time_GNSS(1:end-1))*24*60;
avg_time_res = mean(time_diff);

% time difference between GPS measurements in minutes
gps_times = time_GNSS(sat_index.gps);
time_diff_gps = (gps_times(2:end) - gps_times(1:end-1))*24*60;
avg_time_res_gps = mean(time_diff_gps);

% time difference between GPS+GLONASS measurements in minutes
gps_glo_times = time_GNSS(sat_index.gps_glo);
time_diff_gps_glo = (gps_glo_times(2:end) - gps_glo_times(1:end-1))*24*60;
avg_time_res_gps_glo = mean(time_diff_gps_glo);

%% frequency offset calculation

% load offset 
load('N001_evaluations.mat', 'freq_offset');

% apply frequency offset for each constellation
SL_GNSSIR(sat_index.glonass) = SL_GNSSIR(sat_index.glonass) + freq_offset.glo;

% frequency offset evaluation
[evals2, evaluations2] = compute_Evaluations(num_constellations, SL_GNSSIR_MED, SL_GNSSIR_MEAN, SL_stdev, TG_daily_mean, SL_GNSSIR, SL_TG, sat_index);

% save evaluations
%save([GNSS_station '_evaluations.mat'], 'evals2', 'evaluations2', 'freq_offset', 'SL_GNSSIR','-append');

%% figure 1

%close all;

% GNSS-IR and TG Mean Sea Level Time Series
figure(1);
set(gcf,  'Position', [400, 400, 1300, 300])
plot(time_SL, SL_GNSSIR, 'Linewidth', 1);
hold on; 
grid on;
plot(time_SL, SL_TG, 'Linewidth', 1);
yline(MSL, 'Linewidth', 1.2, 'Color', [227/255,118/255,0/255]);
yline(0, 'Linewidth', 1.2);
xlim([time_SL(1) - 3, time_SL(end) + 3,]);
datetick('x', 'mm/dd', 'keeplimits');
ylabel('Mean Sea Level (m)');
legend('GNSS-IR', 'TG', 'MSL', ['MSL ' TG_epoch], 'Location', 'Southoutside', 'Orientation', 'Horizontal');
title('GNSS-IR vs TG Mean Sea Level Time Series');
subtitle(append('GNSS Station ', GNSS_station, ' and NOAA Station ', TG_station, ', ', datestr(time_GNSS(1), 'mmm yyyy'), ' to ', datestr(time_GNSS(end), 'mmm yyyy')));
set(gca,'fontname','times','FontSize',14);

%% figure 2

%close all; 

% daily mean TGA RH and daily mean and median GNSS-IR RH with stdev
figure(2);
set(gcf,  'Position', [200, 400, 1200, 400])
tgmean_plt = plot(day_GNSS, TG_daily_mean, 'm+', 'Linewidth', 2);
hold on;
grid on;
plot(day_GNSS, TG_daily_mean,'m-', 'Linewidth', 0.5);
xlim([day_GNSS(1) - 3, day_GNSS(end) + 3,]);
med_plt = plot(day_GNSS, SL_GNSSIR_MED,'r*', 'Linewidth', 1);
plot(day_GNSS, SL_GNSSIR_MED,'r-', 'Linewidth', 0.5);
plot(day_GNSS, SL_GNSSIR_MEAN,'b-', 'Linewidth', 0.5);
mean_plt = scatter(day_GNSS, SL_GNSSIR_MEAN, 30, stdev_RH,'filled');
c = colorbar;
datetick('x', 'mm/dd', 'keeplimits');
legend([tgmean_plt, med_plt, mean_plt], {'Tide Gauge Daily Mean', 'GNSS-IR Daily Median', 'GNSS-IR Daily Average'}, 'Location', 'Southeast'); 
title('Tide Gauge Daily Average and GNSS Daily Average and Median Sea Level with Standard Deviations');
subtitle(append('GNSS Station ', GNSS_station, ' and NOAA Station ', TG_station, ', ', datestr(day_GNSS(1), 'mmm yyyy'), ' to ', datestr(day_GNSS(end), 'mmm yyyy')));
ylabel('Mean Sea Level (m)')
ylabel(c, 'Standard Deviation of Daily Average MSL (m)','fontname','times','FontSize',14);
set(gca,'fontname','times','FontSize',14);

%% figure 3

%close all;

start_time = time_SL(round(length(time_SL)/2));
end_time = start_time + 12; 
timelim = [datenum(start_time) datenum(end_time)];

% small section of TG time series with satellite measurements + stdev
figure(3);
set(gcf,  'Position', [400, 400, 1100, 500])
tga_plt = plot(time_TG, TG.Verified_m_, 'Color', [227/255,118/255,0/255], 'Linewidth', 0.75);
hold on;
grid on;
err_plt = errorbar(time_SL, SL_GNSSIR, SL_stdev, '.');
gps_plt = scatter(time_SL(sat_index.gps), SL_GNSSIR(sat_index.gps), 50, satColors(1,:), 'filled');
glo_plt = scatter(time_SL(sat_index.glonass), SL_GNSSIR(sat_index.glonass), 50, satColors(2,:), 'filled');
xlim(timelim);
datetick('x', 'mm/dd', 'keeplimits');
ylabel('Sea Level (m)');
legend([tga_plt, gps_plt, glo_plt, err_plt], {'Tide Gauge', 'GPS', 'GLONASS', 'Uncertainty'}, 'Location', 'Southoutside', 'Orientation', 'Horizontal');
title('Tide Gauge Sea Level Time Series with GNSS-IR Satellite Measurements');
subtitle(append('GNSS Station ', GNSS_station, ' and NOAA Station ', TG_station, ', ', datestr(start_time, 'mmm yyyy')));
set(gca,'fontname','times','FontSize',14);

%% figure 4

%close all;

% pie chart showing distribution of data among satellites
fig = figure(4);
labels = {['GPS (' num2str(round(length(sat_index.gps)/length(RH_GNSS)*100)) '%)'];
    ['GLONASS (' num2str(round(length(sat_index.glonass)/length(RH_GNSS)*100)) '%)'];
    };
pie([length(sat_index.gps), length(sat_index.glonass)], labels);
title(append('Satellite Data Distribution, ', 'GNSS Station ', GNSS_station, ', ', datestr(time_GNSS(1), 'mmm yyyy'), ' to ', datestr(time_GNSS(end), 'mmm yyyy')));
set(gca,'fontname','times','FontSize',14);
set(findobj(fig,'type','text'),'FontName','times','FontSize',14);
ax = gca();
ax.Colormap = satColors(1:num_constellations,:);

%% figure 5

%close all;

% evaluations ranked from lowest to highest RMSD
figure(5);
set(gcf,  'Position', [200, 400, 900, 600])
subplot(2,1,1);
eval1c = plot(1:length(evaluations(:,1)), evaluations(:,1), '-', 'Linewidth', 2, 'Color', '#00a0a2');
hold on;
grid on;
eval2c = plot(1:length(evaluations2(:,1)), evaluations2(:,1), 'x-', 'Linewidth', 2, 'Color', '#00a0a2');
eval1c.Color(4) = 0.4;
ylabel('Coefficient');
title('Correlation');
xlim([0,length(all_labels)+1]);
xticks(1:length(all_labels));
xticklabels(all_labels);
xtickangle(45);
%sgtitle('Evaluation of GNSS-IR Method with Tide Gauge Sea Level', 'fontname','Arial','FontSize', 30, 'Fontweight', 'bold');
legend('Without Offset', 'With Offset');
set(gca,'fontname','Times','FontSize',14);
subplot(2,1,2);
eval1r = plot(1:length(evaluations(:,1)), evaluations(:,2), '-', 'Linewidth', 2, 'Color', '#d3007c');
hold on;
eval2r = plot(1:length(evaluations2(:,1)), evaluations2(:,2), 'x-', 'Linewidth', 2, 'Color', '#d3007c');
eval1r.Color(4) = 0.4;
title('Root-Mean-Square Deviation');
ylabel('RMSD (m)');
grid on;
xlim([0,length(all_labels)+1]);
xticks(1:length(all_labels));
xticklabels(all_labels);
xtickangle(45);
legend('Without Offset', 'With Offset');
set(gca,'fontname','Times','FontSize',14);
sgtitle('Evaluation of GNSS-IR Measurements with an Interpolated Tide Gauge Dataset', 'fontname','Times','FontSize',16);

