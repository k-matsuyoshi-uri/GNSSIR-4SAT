% view variables and QC from SNR data
% Kayli Matsuyoshi

%% setup 
close all;
clear;

%% set parameters
SNR_directory = 'N001_SNR/';
max_azim = 145.1;
min_azim = 15;
max_elev = 16;
min_elev = 5;

%% loop on or off

snr_meta = struct2table(dir(SNR_directory));
files = (snr_meta(3:end,1));

% iterate through files, creating and saving plots and stats
%for f = 1:size(files, 1) % LOOP ON
    
% OR specify file with
f = 1; % LOOP OFF

%% load data 
% set up array of file names used

filename = string(table2cell(files(f,1)));
%filename = 'getting_started/gnss_snr_outputs/npt2019_0329.snr99';
snr = load(append(SNR_directory,filename));

%% level 1 data extraction

% satellite numbers
sat_num = snr(:,1); 

% elevation angle (º)
elevation_angle = snr(:,2);

% azimuth angle (º)
azimuth_angle = snr(:,3);

% time in seconds (GPS day)
time_GPS = snr(:,4);

% elevation angle rate of change (º/s)
elevation_change = snr(:,5);

% S1 (dBHz)
S1 = snr(:,7);

% S2 (dBHz)
S2 = snr(:,8);

%% level 2 data extraction

gps_index = find(sat_num < 100);
glonass_index = find(100 <= sat_num & sat_num < 200);
galileo_index = find(200 <= sat_num & sat_num < 300);
beidou_index = find(300 <= sat_num);

% sort unique satellite numbers from SNR data
all_sats = unique(sat_num);

% identify gps, glonass, galileo, and beidou satellites
gps_sat_num = all_sats(all_sats < 100);
glonass_num = all_sats(100 <= all_sats & all_sats < 200);
galileo_num = all_sats(200 <= all_sats & all_sats < 300);
beidou_num = all_sats(300 <= all_sats);

% elevation angle statistics
[min_elevation, min_elevation_pos] = min(elevation_angle);
[max_elevation, max_elevation_pos] = max(elevation_angle);
mean_elevation = mean(elevation_angle);

% azimuth angle statistics
[min_azimuth, min_azimuth_pos] = min(azimuth_angle);
[max_azimuth, max_azimuth_pos] = max(azimuth_angle);
mean_azimuth = mean(azimuth_angle);

S2_index = find(S2 > 0);
S2_sats = unique(snr(S2_index,1));

% quality control 
azimuth_qc = find(azimuth_angle > min_azim & azimuth_angle < max_azim);
elevation_qc = find(elevation_angle > min_elev & elevation_angle < max_elev);
qc = intersect(azimuth_qc, elevation_qc);

%% plotting data

% no quality control
fig = figure(1);
subplot(1,2,1);
set(gcf, 'Position', [400, 400, 1000, 400]);
scatter(elevation_angle, S1, 10, sat_num, 'filled');
hold on;
grid on;
xlim([5, 30]);
ylim([0, 60]);
xlabel('Elevation Angle (degrees)');
ylabel('SNR (dBHz)');
date = strrep(erase(eraseBetween(filename, 1, 5), '.snr99'), '_', '/');
date = append(extractBefore(date, 8), '/', extractAfter(date, 7));
title(append(extractBefore(filename, 5), ' ', date, ' SNR vs Elevation'));
set(gca,'fontname','times','FontSize',14);
subplot(1,2,2);
scatter(elevation_angle(qc), S1(qc), 10, sat_num(qc), 'filled');
hold on;
grid on;
c = colorbar;
xlim([5, 30]);
ylim([0, 60]);
xlabel('Elevation Angle (degrees)');
ylabel('SNR (dBHz)');
c.Ticks = [20, 120, 220, 320];
labels = {'GPS', 'GLONASS', 'Galileo', 'BeiDou'};
c.TickLabels = labels;
%ylabel(c, 'Satellite Numbers', 'fontname','times','FontSize',14);
title(append(extractBefore(filename, 5), ' ', date, ' QC SNR vs Elevation'));
set(gca,'fontname','times','FontSize',14);
%saveas(fig, append('figures/SNR_vs_elevation/SNR_vs_elevation_QC_',erase(erase(filename, SNR_directory), '.snr99')), 'png');

%% saving data

% save_file = append('SNR_stats/',erase(erase(filename, SNR_directory), '.snr99'),'_SNR_stats.mat');
% %load(save_file);
% used = size(qc)/size(sat_num);
% total_GPS = length(all_sats(all_sats < 100));
% total_GLONASS = length(all_sats(all_sats >= 100 & all_sats < 200));
% total_Galileo = length(all_sats(all_sats >= 200 & all_sats < 300));
% total_Beidou = length(all_sats(all_sats >= 300));
% stats = struct('used', used, 'total_GPS', total_GPS, 'total_GLONASS', total_GLONASS, 'total_Galileo', total_Galileo, 'total_Beidou', total_Beidou, 'max_azimuth', max_azimuth, 'max_azimuth_pos', max_azimuth_pos, 'min_azimuth', min_azimuth, 'min_azimuth_pos', min_azimuth_pos, 'max_elevation', max_elevation, 'max_elevation_pos', max_elevation_pos, 'min_elevation', min_elevation, 'min_elevation_pos', min_elevation_pos);
% SNR = struct('all_sats', all_sats, 'sat_num', sat_num, 'azimuth', azimuth_angle, 'elevation', elevation_angle, 'S1', S1, 'S2', S2, 'S2_sats', S2_sats, 'qc', qc, 'time_GPS', time_GPS);
% save(save_file, 'SNR', 'stats'); % if creating new files

%% close and reset

%close all; % LOOP ON

%end % LOOP ON
%}